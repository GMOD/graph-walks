import {
  RegionTooLargeError,
  bpString,
  panSNHaplotype,
  panSNMatchesPrefix,
} from './util.ts'
import { placeAndLink } from './walkChunk.ts'

import type { GraphTables, WalkStopTable } from './util.ts'
import type { Chunk } from './walkChunk.ts'

// A walk row of gfa-to-tabix's format 5, its steps decoded with its chunk
// (walkChunk.ts). A haplotype path is cut into pieces, one per chunk it
// passes, and each piece's row takes the next piece index, so consecutive
// indices of one path join. A row names the chunk its path's piece before and
// after it is filed under, `prev` and `next`, with none at either end of the
// path. A record's first and last rows also name the gap to the contig's
// record before and after, `prevRecord` and `nextRecord`.
export interface WalkRow {
  name: string
  fragStart: number
  hapOffset: number
  piece: number
  n: number
  // of the n steps, those off R, in routes' interiors
  interior?: number
  steps: { ids: Int32Array; rev: Uint8Array }
  chunkStart: number
  prev?: WalkChunk
  next?: WalkChunk
  prevRecord?: WalkGap
  nextRecord?: WalkGap
  // an `abridged` row, standing for its path's pieces in the chunk: `steps`
  // holds the handles on R it starts and ends on
  abridged?: Abridged
}

/**
 * A path's pieces in a chunk that the writer moved to a detail record: their
 * piece numbers, steps and bp, and the record's sequence and position, which a
 * cut that expands the path reads
 */
export interface Abridged {
  pieces: number[]
  steps: number
  bp: number
  detail: { refName: string; start: number }
}

export interface WalkChunk {
  refName: string
  start: number
}

// Contig bp the graph lacks between two records, and the chunk the other
// record's nearest piece is filed under, when that record has rows
export interface WalkGap {
  bp: number
  chunk?: WalkChunk
}

/**
 * What lies past one end of a fragment. `trimmed` bp of the walk come first
 * where the cut read them but kept none, then:
 * - `elsewhere`: the walk goes on in `chunk`, which the cut did not read; `bp`
 *   is how far on it comes back, when the read holds its next piece
 * - `break`: the contig's record ends, and its next record starts `bp` on;
 *   `chunk` files that record's nearest piece, when it has rows
 * - `end`: the contig's last record in the graph ends
 */
export type WalkStop = (
  | { kind: 'elsewhere'; chunk: WalkChunk; bp?: number }
  | { kind: 'break'; bp: number; chunk?: WalkChunk }
  | { kind: 'end' }
) & { trimmed?: number }

const WALK_STOP_KINDS = ['elsewhere', 'break', 'end'] as const

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

// A pass of the reference over a node longer than a chunk, from a `#long` line
export interface LongPass {
  refName: string
  start: number
  end: number
}

export interface WalkHeader {
  // format:i:5; the reader refuses any other, and format 1 states none
  format?: number
  chunk?: number
  maxNode?: number
  references: string[]
  haplotypes: string[]
  long: LongPass[]
}

function integerTag(line: string, tag: string) {
  const match = new RegExp(`(?:^|\\t)${tag}:i:(\\d+)`).exec(line)
  return match ? +match[1]! : undefined
}

/**
 * What a walk file's header states: a `#walks` line with the format, the
 * chunk size and the longest node in bp, `format:i:5\tchunk:i:65536\tmaxnode:i:1024`,
 * a `#reference\tGRCh38` line naming the reference sample and a
 * `#haplotype\tHG002#1` line per other haplotype with rows
 */
export function walkHeader(lines: string[]): WalkHeader {
  const header: WalkHeader = { references: [], haplotypes: [], long: [] }
  for (const line of lines) {
    const [key, value = '', start = '', end = ''] = line.split('\t')
    if (key === '#reference') {
      header.references.push(value)
    } else if (key === '#haplotype') {
      header.haplotypes.push(value)
    } else if (key === '#long') {
      header.long.push({ refName: value, start: +start, end: +end })
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

/**
 * The first base a cut of `refName` over `start`–`end` queries: the chunk
 * before the window's first, or the chunk before a `#long` pass overlapping
 * the window, whichever is earlier. A file with no `#long` lines reads
 * lookbackChunks back, which is one chunk unless `maxnode` exceeds a chunk,
 * as in a pggb file from gfa-to-tabix 1.2.
 */
export function walkQueryStart(
  header: Pick<WalkHeader, 'maxNode' | 'long'> & { chunk: number },
  refName: string,
  start: number,
  end: number,
) {
  const { chunk, long } = header
  if (long.length === 0) {
    return chunkQueryStart(start, chunk, lookbackChunks(header.maxNode, chunk))
  }
  let from = chunkQueryStart(start, chunk)
  for (const pass of long) {
    if (pass.refName === refName && pass.start < end && pass.end > start) {
      from = Math.min(from, chunkQueryStart(pass.start, chunk))
    }
  }
  return from
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

/**
 * The walks a cut on `refName` decodes: its own and those `wanted` names. A
 * pggb P-line names a reference with no haplotype field (`chm13#chr1`), so the
 * filter matches the reference by its whole name.
 */
export function walksKept(wanted: string[] | undefined, refName: string) {
  return walkNameFilter(wanted, panSNHaplotype(refName) ?? refName)
}

// One stretch of a haplotype path: rows of consecutive piece index, decoded,
// and where the path goes on before and after it, when the rows say
export interface WalkRun {
  name: string
  fragStart: number
  hapOffset: number
  firstPiece: number
  ids: Int32Array
  rev: Uint8Array
  prev?: WalkChunk
  next?: WalkChunk
  prevRecord?: WalkGap
  nextRecord?: WalkGap
  abridged?: Abridged
}

export function joinPieces(rows: WalkRow[]) {
  const byPath = new Map<string, WalkRow[]>()
  const runs: WalkRun[] = []
  for (const row of rows) {
    if (row.abridged) {
      runs.push(decodeRun([row]))
      continue
    }
    const key = `${row.name}\t${row.fragStart}`
    const list = byPath.get(key)
    if (list) {
      list.push(row)
    } else {
      byPath.set(key, [row])
    }
  }
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
  const last = pieces.at(-1)!
  return {
    name: first.name,
    fragStart: first.fragStart,
    hapOffset: first.hapOffset,
    firstPiece: first.piece,
    ids,
    rev,
    prev: first.prev,
    next: last.next,
    prevRecord: first.prevRecord,
    nextRecord: last.nextRecord,
    ...(first.abridged ? { abridged: first.abridged } : {}),
  }
}

export interface WalkFragment {
  name: string
  fragStart: number
  hapStart: number
  hapEnd: number
  ids: Int32Array
  rev: Uint8Array
  before: WalkStop
  after: WalkStop
  // set on the fragment of an abridged row, from where it leaves R to where it
  // rejoins, which runs `abridged.bp` on the haplotype
  abridged?: Abridged
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
 *
 * A contained cut keeps the reference over the window and its context alone,
 * and splits a run wherever it steps onto the reference off that stretch.
 * Each piece between those steps that steps on the stretch is kept from its
 * first such step to its last, so a detour that comes back is kept whole, and
 * `context` bp of the walk past either end. At AMY1 haplotypes go from a 4 kb
 * window to the other copies and back: the cut above keeps 9,928 nodes there,
 * the contained cut 182.
 *
 * Each fragment says what lies past either end (WalkStop).
 */
export function walkCut(
  runs: WalkRun[],
  nodes: Map<number, WalkNode>,
  refName: string,
  { contained = false, context = 0 } = {},
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
  const length = (id: number) => {
    const node = nodes.get(id)
    return node === undefined ? 0 : node.end - node.start
  }
  // A contained run's spans: each stretch between steps onto the reference
  // off the window that has steps on it, from its first such step to its
  // last, and `context` bp of the walk past either
  const containedSpans = (run: WalkRun) => {
    const { ids } = run
    let from = 0
    let first = -1
    let last = -1
    for (let i = 0; i <= ids.length; i++) {
      const node = i < ids.length ? nodes.get(ids[i]!) : undefined
      if (i === ids.length || (node?.rank === 0 && !node.onWindow)) {
        if (first >= 0) {
          for (let bp = 0; first > from && bp < context;) {
            bp += length(ids[--first]!)
          }
          for (let bp = 0; last < i - 1 && bp < context;) {
            bp += length(ids[++last]!)
          }
          if (first === 0) {
            goesOn(run.prev)
          }
          if (last === ids.length - 1) {
            goesOn(run.next)
          }
          spans.push({ run, first, last })
        }
        from = i + 1
        first = -1
      } else if (node?.onWindow) {
        if (first < 0) {
          first = i
        }
        last = i
      }
    }
  }
  const detourSpan = (run: WalkRun) => {
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
  // the window and its context, as the reference nodes on it span them
  let windowStart = Infinity
  let windowEnd = -Infinity
  for (const node of nodes.values()) {
    if (node.onWindow) {
      windowStart = Math.min(windowStart, node.start)
      windowEnd = Math.max(windowEnd, node.end)
    }
  }
  // an abridged run is kept where its stretch of R overlaps the window
  const abridgedSpan = (run: WalkRun) => {
    const a = reference(run.ids[0]!)
    const b = reference(run.ids.at(-1)!)
    if (a?.refName !== refName || b?.refName !== refName) {
      return
    }
    if (
      Math.max(a.end, b.end) > windowStart &&
      Math.min(a.start, b.start) < windowEnd
    ) {
      spans.push({ run, first: 0, last: run.ids.length - 1 })
    }
  }
  for (const run of runs) {
    if (run.abridged) {
      abridgedSpan(run)
    } else if (contained) {
      containedSpans(run)
    } else {
      detourSpan(run)
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
  const abridgedRuns = new Set(
    spans.filter(span => span.run.abridged).map(span => span.run),
  )
  const bpOf = (ids: Int32Array, from: number, to: number) => {
    let bp = 0
    for (let i = from; i < to; i++) {
      bp += length(ids[i]!)
    }
    return bp
  }
  // each run's neighbours in its record, by piece index
  const byRecord = new Map<string, WalkRun[]>()
  for (const run of runs) {
    const key = `${run.name}\t${run.fragStart}`
    const list = byRecord.get(key)
    if (list) {
      list.push(run)
    } else {
      byRecord.set(key, [run])
    }
  }
  const neighbours = new Map<WalkRun, [WalkRun?, WalkRun?]>()
  for (const list of byRecord.values()) {
    list.sort((a, b) => a.firstPiece - b.firstPiece)
    list.forEach((run, i) => {
      neighbours.set(run, [list[i - 1], list[i + 1]])
    })
  }
  const ends = new Map<WalkRun, number>()
  const runEnd = (run: WalkRun) => {
    let end = ends.get(run)
    if (end === undefined) {
      end =
        run.hapOffset + (run.abridged?.bp ?? bpOf(run.ids, 0, run.ids.length))
      ends.set(run, end)
    }
    return end
  }
  const stop = (run: WalkRun, after: boolean, trimmed: number): WalkStop => {
    const cut = trimmed > 0 ? { trimmed } : {}
    const chunk = after ? run.next : run.prev
    if (chunk) {
      const other = neighbours.get(run)?.[after ? 1 : 0]
      return other
        ? {
            kind: 'elsewhere',
            chunk,
            bp: after
              ? other.hapOffset - runEnd(run)
              : run.hapOffset - runEnd(other),
            ...cut,
          }
        : { kind: 'elsewhere', chunk, ...cut }
    }
    const gap = after ? run.nextRecord : run.prevRecord
    return gap ? { kind: 'break', ...gap, ...cut } : { kind: 'end', ...cut }
  }
  const fragments: WalkFragment[] = []
  for (const { run, first: f, last: l } of spans) {
    const { ids, rev } = run
    if (abridgedRuns.has(run)) {
      fragments.push({
        name: run.name,
        fragStart: run.fragStart,
        hapStart: run.hapOffset,
        hapEnd: runEnd(run),
        ids,
        rev,
        before: stop(run, false, 0),
        after: stop(run, true, 0),
        abridged: run.abridged,
      })
      continue
    }
    let first = f
    let last = l
    while (first > 0 && kept.has(ids[first - 1]!)) {
      first--
    }
    while (last < ids.length - 1 && kept.has(ids[last + 1]!)) {
      last++
    }
    const hapStart = run.hapOffset + bpOf(ids, 0, first)
    const hapEnd = hapStart + bpOf(ids, first, last + 1)
    fragments.push({
      name: run.name,
      fragStart: run.fragStart,
      hapStart,
      hapEnd,
      ids: ids.subarray(first, last + 1),
      rev: rev.subarray(first, last + 1),
      before: stop(run, false, hapStart - run.hapOffset),
      after: stop(run, true, runEnd(run) - hapEnd),
    })
  }
  return { kept, fragments, open: [...open.values()] }
}

/**
 * Adds to `into` the nodes of decoded rows that walkCut keeps whatever else
 * the cut reads: those of `chunk.nodes` on the window, and those a row steps
 * on between two of its steps on the window, with no step onto the reference
 * off the window between them in a contained cut. A node placed later at a
 * reference position off the window can split a contained walk, so the count
 * holds for the placements at the time.
 */
export function surelyKept(
  chunk: { nodes: Map<number, unknown>; rows: WalkRow[] },
  nodes: Map<number, WalkNode>,
  contained: boolean,
  into: Set<number>,
) {
  for (const id of chunk.nodes.keys()) {
    if (nodes.get(id)?.onWindow) {
      into.add(id)
    }
  }
  for (const { steps } of chunk.rows) {
    const { ids } = steps
    let last = -1
    for (let i = 0; i < ids.length; i++) {
      const node = nodes.get(ids[i]!)
      if (node?.onWindow) {
        for (let j = last + 1; last >= 0 && j < i; j++) {
          into.add(ids[j]!)
        }
        last = i
      } else if (contained && node?.rank === 0) {
        last = -1
      }
    }
  }
}

function stopTable(stops: WalkStop[]): WalkStopTable {
  const chunkRefNames: string[] = []
  const refIndex = new Map<string, number>()
  const table: WalkStopTable = {
    kinds: new Uint8Array(stops.length),
    kindNames: [...WALK_STOP_KINDS],
    bp: new Float64Array(stops.length).fill(Number.NaN),
    trimmed: new Float64Array(stops.length),
    chunkRefs: new Int32Array(stops.length).fill(-1),
    chunkStarts: new Float64Array(stops.length).fill(Number.NaN),
    chunkRefNames,
  }
  stops.forEach((stop, i) => {
    table.kinds[i] = WALK_STOP_KINDS.indexOf(stop.kind)
    table.trimmed[i] = stop.trimmed ?? 0
    if ('bp' in stop && stop.bp !== undefined) {
      table.bp[i] = stop.bp
    }
    if ('chunk' in stop && stop.chunk) {
      const { refName, start } = stop.chunk
      let ref = refIndex.get(refName)
      if (ref === undefined) {
        ref = chunkRefNames.length
        refIndex.set(refName, ref)
        chunkRefNames.push(refName)
      }
      table.chunkRefs[i] = ref
      table.chunkStarts[i] = start
    }
  })
  return table
}

// Whether placement `node` beats `held`: a reference position beats any
// other, and among those the path sorting first at its earliest start wins,
// so the order a cut reads chunks in does not matter
const outranks = (
  node: Omit<WalkNode, 'onWindow'>,
  held: Omit<WalkNode, 'onWindow'>,
) =>
  node.rank === 0 &&
  (held.rank !== 0 ||
    node.refName < held.refName ||
    (node.refName === held.refName && node.start < held.start))

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

  // A node with its position, unless the cut already holds it at one that
  // outranks it
  placeNode(id: number, node: Omit<WalkNode, 'onWindow'>) {
    const held = this.nodes.get(id)
    if (held && !outranks(node, held)) {
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
        before: stopTable(walks.map(f => f.before)),
        after: stopTable(walks.map(f => f.after)),
        abridgedSteps: Float64Array.from(walks, f => f.abridged?.steps ?? NaN),
        abridgedBp: Float64Array.from(walks, f => f.abridged?.bp ?? NaN),
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

// The reader decodes a step along R several times as fast as one off it, so
// the read budget counts it as a quarter
const ALONG_R_PER_STEP = 4

const readSteps = (row: WalkRow) => {
  const interior = row.interior ?? 0
  return interior + (row.n - interior) / ALONG_R_PER_STEP
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
    const steps = readSteps(row)
    total += steps
    perChunk.set(row.chunkStart, (perChunk.get(row.chunkStart) ?? 0) + steps)
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
      : `Too many haplotype steps here to draw ${filtered ? 'these haplotypes' : 'every haplotype'} (${Math.round(read).toLocaleString()} steps read, a step along the reference counting a quarter, over ${readPerKept} × ${names.steps} of ${stepBudget.toLocaleString()}); choose fewer haplotypes`,
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
  return keptOverError(total, budget, window, filtered)
}

// keptStepsError's notice for a cut that keeps more nodes than `limit`, which
// it counted on the window before `readTo` when it stopped reading there
export function keptNodesError(
  nodes: number,
  limit: number,
  window: { start: number; end: number },
  filtered: boolean,
  readTo = window.end,
) {
  return keptOverError(nodes, limit, window, filtered, readTo)
}

function keptOverError(
  kept: number,
  limit: number,
  window: { start: number; end: number },
  filtered: boolean,
  readTo = window.end,
) {
  if (kept <= limit) {
    return undefined
  }
  const windowBp = window.end - window.start
  const readBp = Math.min(readTo, window.end) - window.start
  const fitsBp = Math.max(1, Math.floor((readBp * limit) / kept))
  return new RegionTooLargeError(
    `Zoom in to about ${bpString(fitsBp)}${filtered ? '' : ', or choose fewer haplotypes,'} to see the graph`,
    limit,
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
  return fitsBp > 0
    ? new RegionTooLargeError(
        `Zoom in to about ${bpString(fitsBp)} to see the graph`,
        budget,
        windowBp,
        fitsBp,
      )
    : new RegionTooLargeError(
        `Too much graph here to fetch, even zoomed in (${megabytes(least ?? bytes)} against ${names.bytes} ${megabytes(budget)}). A satellite or rDNA array usually causes this: the graph folds the repeat into a few nodes that each haplotype loops through many times`,
        budget,
        windowBp,
        1,
        { bytes },
      )
}

const megabytes = (bytes: number) => `${(bytes / 1e6).toFixed(1)} MB`
