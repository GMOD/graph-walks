import { RegionTooLargeError, bpString, panSNMatchesPrefix } from './util.ts'
import { placeAndLink } from './walkChunk.ts'

import type { GraphTables } from './util.ts'
import type { Chunk } from './walkChunk.ts'

// A walk row of gfa-to-tabix's format 2, its steps decoded with its chunk
// (walkFormat2.ts). A haplotype path is cut into pieces, one per chunk it
// passes, and each piece's row takes the next piece index, so consecutive
// indices of one path join. A row names the chunk its path's piece before and
// after it is filed under, `prev` and `next`, with none at either end of the
// path.
export interface WalkRow {
  name: string
  fragStart: number
  hapOffset: number
  piece: number
  n: number
  steps: { ids: Int32Array; rev: Uint8Array }
  chunkStart: number
  prev?: WalkChunk
  next?: WalkChunk
}

export interface WalkChunk {
  refName: string
  start: number
}

export interface WalkNode {
  refName: string
  start: number
  end: number
  rank: number
  // on the reference inside the window plus the cut's context
  onWindow: boolean
}

export interface WalkWindow {
  refName: string
  start: number
  end: number
  context: number
  // keep each node's bases, from the SQ:Z: tags `--sequences` writes
  sequences?: boolean
}

export interface WalkHeader {
  // 2 for gfa-to-tabix's format 2 (format:i:2); format 1 states none, and
  // the reader refuses it
  format?: number
  chunk?: number
  maxNode?: number
  references: string[]
  haplotypes: string[]
}

function integerTag(line: string, tag: string) {
  const match = new RegExp(`(?:^|\\t)${tag}:i:(\\d+)`).exec(line)
  return match ? +match[1]! : undefined
}

/**
 * What a walk file's header states: a `#walks` line with the format, the
 * chunk size and the longest node in bp, `format:i:2\tchunk:i:65536\tmaxnode:i:1024`,
 * a `#reference\tGRCh38` line naming the reference sample and a
 * `#haplotype\tHG002#1` line per other haplotype with rows
 */
export function walkHeader(lines: string[]): WalkHeader {
  const header: WalkHeader = { references: [], haplotypes: [] }
  for (const line of lines) {
    const [key, value = ''] = line.split('\t')
    if (key === '#reference') {
      header.references.push(value)
    } else if (key === '#haplotype') {
      header.haplotypes.push(value)
    } else {
      header.format ??= integerTag(line, 'format')
      header.chunk ??= integerTag(line, 'chunk')
      header.maxNode ??= integerTag(line, 'maxnode')
    }
  }
  return header
}

/**
 * How many chunks before the window's own a cut reads. A node is filed under
 * the chunk holding its start, so a node `maxNode` bp long reaching into the
 * window starts at most ceil(maxNode / chunk) chunks back. At least one, since
 * the chunk before carries the detour a haplotype takes into the window's
 * first reference node; one when the header gives no `maxnode:i:`, which holds
 * for vg's graphs, chopped to nodes of at most 1,024 bp.
 */
export function lookbackChunks(maxNode: number | undefined, chunk: number) {
  return maxNode === undefined ? 1 : Math.max(1, Math.ceil(maxNode / chunk))
}

// The first base a cut queries: the start of the chunk `lookback` chunks
// before the one holding the window's start
export function chunkQueryStart(start: number, chunk: number, lookback = 1) {
  return Math.max(0, (Math.floor(start / chunk) - lookback) * chunk)
}

// Which walks a cut decodes: the reference's, and those `wanted` names by
// PanSN prefix at sample or haplotype depth; undefined decodes every one
export function walkNameFilter(
  wanted: string[] | undefined,
  referenceHaplotype: string | undefined,
) {
  return wanted === undefined || wanted.length === 0
    ? undefined
    : (name: string) =>
        panSNMatchesPrefix(name, referenceHaplotype) ||
        wanted.some(prefix => panSNMatchesPrefix(name, prefix))
}

// One stretch of a haplotype path: rows of consecutive piece index, decoded,
// and where the path goes on before and after it, when the rows say
export interface WalkRun {
  name: string
  hapOffset: number
  ids: Int32Array
  rev: Uint8Array
  prev?: WalkChunk
  next?: WalkChunk
}

export function joinPieces(rows: WalkRow[]) {
  const byPath = new Map<string, WalkRow[]>()
  for (const row of rows) {
    const key = `${row.name}\t${row.fragStart}`
    const list = byPath.get(key)
    if (list) {
      list.push(row)
    } else {
      byPath.set(key, [row])
    }
  }
  const runs: WalkRun[] = []
  for (const pieces of byPath.values()) {
    pieces.sort((a, b) => a.piece - b.piece)
    let from = 0
    for (let i = 1; i <= pieces.length; i++) {
      if (
        i === pieces.length ||
        pieces[i]!.piece !== pieces[i - 1]!.piece + 1
      ) {
        runs.push(decodeRun(pieces.slice(from, i)))
        from = i
      }
    }
  }
  return runs
}

function decodeRun(pieces: WalkRow[]): WalkRun {
  const n = pieces.reduce((sum, p) => sum + p.n, 0)
  const ids = new Int32Array(n)
  const rev = new Uint8Array(n)
  let at = 0
  for (const p of pieces) {
    ids.set(p.steps.ids, at)
    rev.set(p.steps.rev, at)
    at += p.n
  }
  const first = pieces[0]!
  return {
    name: first.name,
    hapOffset: first.hapOffset,
    ids,
    rev,
    prev: first.prev,
    next: pieces.at(-1)!.next,
  }
}

export interface WalkFragment {
  name: string
  hapStart: number
  hapEnd: number
  ids: Int32Array
  rev: Uint8Array
}

/**
 * The cut a window draws from its walks: the reference inside the window plus
 * its context, and each run from its first to its last step on that stretch.
 * A run that leaves the reference there, through steps off it or across a
 * deletion, is followed to the reference step where it rejoins, when the read
 * holds it, and otherwise through the steps off the reference it holds: at
 * the ABCA7 VNTR four haplotypes leave inside the repeat and rejoin 1.8 kb
 * past the window, and cut at the context they measured as partial walks. The reference is cut out as far as those rejoins reach, and
 * a run then runs on outward while the next step's node is already in the
 * cut. Without the run-on a haplotype whose alternate allele straddles the
 * edge stopped short of a node other walks had brought into the cut, and that
 * node drew without it. The run-on adds no node, so one pass is enough, but
 * it reaches through nodes other walks brought in, so with a window edge inside
 * a collapsed array a row's length depends on which haplotypes are cut.
 */
export function walkCut(
  runs: WalkRun[],
  nodes: Map<number, WalkNode>,
  refName: string,
) {
  const onWindow = (id: number) => nodes.get(id)?.onWindow === true
  const reference = (id: number) => {
    const node = nodes.get(id)
    return node?.rank === 0 ? node : undefined
  }
  // the chunks walks cut off at the end of a run go on in
  const open = new Map<string, WalkChunk>()
  const goesOn = (chunk: WalkChunk | undefined) => {
    if (chunk) {
      open.set(`${chunk.refName}:${chunk.start}`, chunk)
    }
  }
  // Where a walk stops being followed from step `at`, going by `dir`: the
  // reference step it rejoins at past any steps off the reference, or across
  // a deletion; else its last step off the reference before the run ends,
  // noting the chunk the path goes on in, or reaches another reference
  // sequence
  const rejoin = (run: WalkRun, at: number, dir: number) => {
    const { ids } = run
    let i = at + dir
    while (i >= 0 && i < ids.length && reference(ids[i]!) === undefined) {
      i += dir
    }
    if (i < 0 || i >= ids.length) {
      goesOn(dir < 0 ? run.prev : run.next)
      return i - dir
    }
    const from = reference(ids[at]!)!
    const to = reference(ids[i]!)!
    if (to.refName !== refName) {
      return i - dir
    }
    const adjacent =
      i === at + dir && (to.start === from.end || to.end === from.start)
    return adjacent ? at : i
  }
  let reachStart = Infinity
  let reachEnd = -Infinity
  const spans: { run: WalkRun; first: number; last: number }[] = []
  for (const run of runs) {
    const { ids } = run
    let first = -1
    let last = -1
    for (let i = 0; i < ids.length; i++) {
      if (onWindow(ids[i]!)) {
        if (first < 0) {
          first = i
        }
        last = i
      }
    }
    if (first >= 0) {
      if (run.name !== refName) {
        first = rejoin(run, first, -1)
        last = rejoin(run, last, 1)
        for (const i of [first, last]) {
          const node = reference(ids[i]!)
          if (node) {
            reachStart = Math.min(reachStart, node.start)
            reachEnd = Math.max(reachEnd, node.end)
          }
        }
      }
      spans.push({ run, first, last })
    }
  }
  const inReach = (id: number) => {
    const node = reference(id)
    return (
      onWindow(id) ||
      (node?.refName === refName &&
        node.end > reachStart &&
        node.start < reachEnd)
    )
  }
  const kept = new Set<number>()
  for (const [id, node] of nodes) {
    if (node.onWindow) {
      kept.add(id)
    }
  }
  for (const span of spans) {
    const { ids } = span.run
    if (span.run.name === refName) {
      while (span.first > 0 && inReach(ids[span.first - 1]!)) {
        span.first--
      }
      while (span.last < ids.length - 1 && inReach(ids[span.last + 1]!)) {
        span.last++
      }
    }
    for (let i = span.first; i <= span.last; i++) {
      kept.add(ids[i]!)
    }
  }
  const length = (id: number) => {
    const node = nodes.get(id)
    return node === undefined ? 0 : node.end - node.start
  }
  const fragments: WalkFragment[] = []
  for (const { run, first: f, last: l } of spans) {
    const { ids, rev } = run
    let first = f
    let last = l
    while (first > 0 && kept.has(ids[first - 1]!)) {
      first--
    }
    while (last < ids.length - 1 && kept.has(ids[last + 1]!)) {
      last++
    }
    let hapStart = run.hapOffset
    for (let i = 0; i < first; i++) {
      hapStart += length(ids[i]!)
    }
    let hapEnd = hapStart
    for (let i = first; i <= last; i++) {
      hapEnd += length(ids[i]!)
    }
    fragments.push({
      name: run.name,
      hapStart,
      hapEnd,
      ids: ids.subarray(first, last + 1),
      rev: rev.subarray(first, last + 1),
    })
  }
  return { kept, fragments, open: [...open.values()] }
}

/**
 * The nodes and links of the chunks a cut reads, and the tables of the cut
 * they make with its walks. A node in two chunks is placed once.
 */
export class WalkGraph {
  readonly nodes = new Map<number, WalkNode>()
  readonly sequences = new Map<number, string>()
  private links: number[] = []
  // the chunks read, whose nodes off R and links tables() fills in
  private readonly chunks: Chunk[] = []

  private readonly window: WalkWindow

  constructor(window: WalkWindow) {
    this.window = window
  }

  // A node with its position, unless the cut already holds it
  placeNode(id: number, node: Omit<WalkNode, 'onWindow'>) {
    if (this.nodes.has(id)) {
      return
    }
    const { window } = this
    this.nodes.set(id, {
      ...node,
      onWindow:
        node.rank === 0 &&
        node.refName === window.refName &&
        node.end > window.start - window.context &&
        node.start < window.end + window.context,
    })
  }

  // A link between two endpoints, each id * 2 + 1 when reversed
  addLinkEnds(source: number, target: number) {
    this.links.push(source, target)
  }

  addChunk(chunk: Chunk) {
    this.chunks.push(chunk)
  }

  // Bases of nodes a chunk's rows give, for a window that asks for them
  addSequences(sequences: Map<number, string>) {
    if (this.window.sequences) {
      for (const [id, bases] of sequences) {
        this.sequences.set(id, bases)
      }
    }
  }

  /**
   * The cut as typed arrays, which the worker transfers rather than copies:
   * nodes by reference position, each link once in the order its L line
   * sorts, then a walk per fragment, the reference's first and the rest by
   * name. graphTablesGFA writes them as the GFA they stand for, the same bytes
   * for the same rows.
   */
  tables(kept: Set<number>, fragments: WalkFragment[]): GraphTables {
    if (this.chunks.length > 0) {
      placeAndLink(
        this.chunks,
        kept,
        (id, node) => {
          const held = this.nodes.get(id)
          if (held) {
            this.nodes.set(id, { ...held, ...node })
          }
        },
        (a, b) => {
          this.links.push(a, b)
        },
      )
    }
    const held: (WalkNode & { id: number })[] = []
    for (const id of kept) {
      const node = this.nodes.get(id)
      if (node) {
        held.push({ ...node, id })
      }
    }
    held.sort((x, y) =>
      x.refName < y.refName
        ? -1
        : x.refName > y.refName
          ? 1
          : x.start - y.start || x.id - y.id,
    )
    const index = new Map<number, number>()
    const names: string[] = []
    const refNames: string[] = []
    const refIndex = new Map<string, number>()
    const lengths = new Int32Array(held.length)
    const refs = new Int32Array(held.length)
    const starts = new Float64Array(held.length)
    const ranks = new Int32Array(held.length)
    held.forEach(({ id, refName, start, end, rank }, i) => {
      index.set(id, i)
      names.push(String(id))
      let ref = refIndex.get(refName)
      if (ref === undefined) {
        ref = refNames.length
        refIndex.set(refName, ref)
        refNames.push(refName)
      }
      lengths[i] = end - start
      refs[i] = ref
      starts[i] = start
      ranks[i] = rank
    })
    // a kept node with no row of its own is named, and not drawn
    const nodeIndex = (id: number) => {
      let i = index.get(id)
      if (i === undefined) {
        i = names.length
        index.set(id, i)
        names.push(String(id))
      }
      return i
    }

    const byLine = new Map<string, number>()
    for (let i = 0; i < this.links.length; i += 2) {
      const a = this.links[i]!
      const b = this.links[i + 1]!
      const source = Math.floor(a / 2)
      const target = Math.floor(b / 2)
      if (kept.has(source) && kept.has(target)) {
        byLine.set(
          `${source}\t${a % 2 ? '-' : '+'}\t${target}\t${b % 2 ? '-' : '+'}`,
          i,
        )
      }
    }
    const lines = [...byLine.keys()].sort()
    const from = new Int32Array(lines.length)
    const to = new Int32Array(lines.length)
    const strands = new Uint8Array(lines.length)
    lines.forEach((line, k) => {
      const i = byLine.get(line)!
      const a = this.links[i]!
      const b = this.links[i + 1]!
      from[k] = nodeIndex(Math.floor(a / 2))
      to[k] = nodeIndex(Math.floor(b / 2))
      strands[k] = (a % 2) | ((b % 2) << 1)
    })

    const walks = this.referenceFirst(fragments)
    const offsets = new Int32Array(walks.length + 1)
    walks.forEach(({ ids }, p) => {
      offsets[p + 1] = offsets[p]! + ids.length
    })
    const steps = new Int32Array(offsets[walks.length]!)
    const reversed = new Uint8Array(steps.length)
    walks.forEach(({ ids, rev }, p) => {
      const at = offsets[p]!
      for (let j = 0; j < ids.length; j++) {
        steps[at + j] = nodeIndex(ids[j]!)
      }
      reversed.set(rev, at)
    })
    return {
      nodes: { names, lengths, refs, starts, ranks, refNames },
      links: { from, to, strands },
      walks: {
        names: walks.map(f => f.name),
        starts: Float64Array.from(walks, f => f.hapStart),
        ends: Float64Array.from(walks, f => f.hapEnd),
        offsets,
        steps,
        reversed,
      },
    }
  }

  // The view anchors on the first walk when no walk names the assembly, as
  // gbz-base's cut has it; in file order that was whichever haplotype the
  // chunk listed first, CHM13 on a GRCh38 window
  private referenceFirst(fragments: WalkFragment[]) {
    const { refName } = this.window
    return [...fragments].sort(
      (a, b) =>
        Number(b.name === refName) - Number(a.name === refName) ||
        (a.name < b.name ? -1 : a.name > b.name ? 1 : a.hapStart - b.hapStart),
    )
  }
}

/**
 * The zoom-in notice for a cut whose walk rows hold more steps than
 * `readPerKept` times `stepBudget`, or undefined when they fit. A window reads whole chunks, those before it
 * included, so the span that fits is counted in chunks from the window's
 * start; when the chunks a window there always reads are over on their own,
 * no zoom fits and the notice asks for fewer haplotypes. A cut for every
 * haplotype is also offered fewer where the reference and one haplotype, at
 * the walks' average steps, would fit the whole window.
 */
// How a notice names the budgets, such as the config slots that set them
export interface BudgetNames {
  steps: string
  bytes: string
}

export const BUDGET_NAMES: BudgetNames = {
  steps: 'the step budget',
  bytes: 'a byte budget of',
}

export function stepBudgetError(
  rows: WalkRow[],
  stepBudget: number,
  readPerKept: number,
  window: { start: number; end: number },
  chunk: number,
  filtered: boolean,
  names = BUDGET_NAMES,
) {
  let total = 0
  const perChunk = new Map<number, number>()
  const walks = new Set<string>()
  for (const row of rows) {
    total += row.n
    perChunk.set(row.chunkStart, (perChunk.get(row.chunkStart) ?? 0) + row.n)
    walks.add(row.name)
  }
  const budget = stepBudget * readPerKept
  if (total <= budget) {
    return undefined
  }
  const windowBp = window.end - window.start
  let read = 0
  let fitsBp = 0
  for (const cs of [...perChunk.keys()].sort((a, b) => a - b)) {
    read += perChunk.get(cs)!
    if (read > budget) {
      break
    }
    fitsBp = Math.max(fitsBp, cs + chunk - window.start)
  }
  const fewerFit =
    !filtered && walks.size > 2 && (2 * total) / walks.size <= budget
  return new RegionTooLargeError(
    fitsBp > 0
      ? `Zoom in to about ${bpString(fitsBp)}${fewerFit ? ', or choose fewer haplotypes,' : ''} to see the graph`
      : `Too many haplotype steps here to draw ${filtered ? 'these haplotypes' : 'every haplotype'} (${read.toLocaleString()} steps read, over ${readPerKept} × ${names.steps} of ${stepBudget.toLocaleString()}); choose fewer haplotypes`,
    budget,
    windowBp,
    Math.max(fitsBp, 1),
  )
}

/**
 * The zoom-in notice for a cut that keeps more steps than `budget`, or
 * undefined when it fits. The span that fits is the window scaled by how far
 * over the cut is.
 */
export function keptStepsError(
  fragments: WalkFragment[],
  budget: number,
  window: { start: number; end: number },
  filtered: boolean,
) {
  const total = fragments.reduce((sum, f) => sum + f.ids.length, 0)
  if (total <= budget) {
    return undefined
  }
  const windowBp = window.end - window.start
  const fitsBp = Math.max(1, Math.floor((windowBp * budget) / total))
  return new RegionTooLargeError(
    `Zoom in to about ${bpString(fitsBp)}${filtered ? '' : ', or choose fewer haplotypes,'} to see the graph`,
    budget,
    windowBp,
    fitsBp,
  )
}

/**
 * The zoom-in notice for a cut whose three files would fetch more than
 * `budget` compressed bytes, or undefined when they fit. `bytesTo(end)` is the
 * indexes' estimate for a read from the window's first query base to `end`;
 * the span that fits grows a chunk at a time from the chunk holding the
 * window's start, and when that chunk and those before it are over on their
 * own, no zoom fits.
 */
export async function byteBudgetError(
  bytesTo: (end: number) => Promise<number>,
  budget: number,
  window: { start: number; end: number },
  chunk: number,
  names = BUDGET_NAMES,
) {
  const bytes = await bytesTo(window.end)
  if (bytes <= budget) {
    return undefined
  }
  let fitsBp = 0
  let least: number | undefined
  for (
    let cs = Math.floor(window.start / chunk) * chunk;
    cs < window.end;
    cs += chunk
  ) {
    const read = await bytesTo(cs + 1)
    least ??= read
    if (read > budget) {
      break
    }
    fitsBp = cs + chunk - window.start
  }
  const windowBp = window.end - window.start
  return new RegionTooLargeError(
    fitsBp > 0
      ? `Zoom in to about ${bpString(fitsBp)} to see the graph`
      : `Too much graph here to fetch (${megabytes(least ?? bytes)} against ${names.bytes} ${megabytes(budget)})`,
    budget,
    windowBp,
    Math.max(fitsBp, 1),
  )
}

const megabytes = (bytes: number) => `${(bytes / 1e6).toFixed(1)} MB`
