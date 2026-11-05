import {
  RegionTooLargeError,
  bpString,
  panSNHaplotype,
  panSNMatchesPrefix,
} from './util.ts'
import { placeAndLink } from './walkChunk.ts'
import { ATTACH_BP, departureStretch } from './walkDeparture.ts'

import type { DepartureTable, GraphTables } from './util.ts'
import type { Chunk } from './walkChunk.ts'
import type { Hole, WalkDeparture } from './walkDeparture.ts'

// A walk row of gfa-to-tabix's format 6, its steps decoded with its chunk
// (walkChunk.ts). A haplotype path is cut into pieces, one per chunk it
// passes, and each piece's row takes the next piece index, so consecutive
// indices of one path join.
export interface WalkRow {
  name: string
  fragStart: number
  hapOffset: number
  piece: number
  n: number
  // of the n steps, those off R, in routes' interiors
  interior?: number
  // refStarts: where each step lies on the reference, NaN off R or unknown;
  // against: 1 where the step runs against R's orientation at its pass
  steps: {
    ids: Int32Array
    rev: Uint8Array
    refStarts?: Float64Array
    against?: Uint8Array
  }
  chunkStart: number
  // an `abridged` row, standing for one piece: `steps` holds the handles on
  // R it starts and ends on
  abridged?: Abridged
}

/**
 * A piece the writer moved to a detail record: all its steps, which reading
 * the record decodes; the bp of its stretch from its first step on R to its
 * last, which runs from its row's hapOffset to hapOffset + bp on the contig;
 * and the record's sequence and position, which a cut that expands the path
 * reads
 */
export interface Abridged {
  steps: number
  bp: number
  detail: { refName: string; start: number }
}

export interface WalkNode {
  refName: string
  start: number
  end: number
  rank: number
  // how many times the reference's paths pass the node, from the ref rows the
  // cut read; 0 for a node it places off R, whether off the reference or on
  // it only where no ref row it read holds the node
  passes: number
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
  // format:i:6; the reader refuses any other, and format 1 states none
  format?: number
  chunk?: number
  maxNode?: number
  references: string[]
  haplotypes: string[]
  long: LongPass[]
  holes: Hole[]
  // how many nodes the reference's paths pass more than once: 0 states
  // every position and departure in the file as a fact of the graph
  revisited?: number
}

function integerTag(line: string, tag: string) {
  const match = new RegExp(`(?:^|\\t)${tag}:i:(\\d+)`).exec(line)
  return match ? +match[1]! : undefined
}

/**
 * What a walk file's header states: a `#walks` line with the format, the
 * chunk size and the longest node in bp, `format:i:6\tchunk:i:65536\tmaxnode:i:1024`,
 * a `#reference\tGRCh38` line naming the reference sample, a
 * `#haplotype\tHG002#1` line per other haplotype with rows, a `#long` line
 * per pass of the reference over a node longer than a chunk and a `#hole`
 * line per stretch of a reference sequence no line of the reference covers
 */
export function walkHeader(lines: string[]): WalkHeader {
  const header: WalkHeader = {
    references: [],
    haplotypes: [],
    long: [],
    holes: [],
  }
  for (const line of lines) {
    const [key, value = '', start = '', end = ''] = line.split('\t')
    if (key === '#reference') {
      header.references.push(value)
    } else if (key === '#haplotype') {
      header.haplotypes.push(value)
    } else if (key === '#long') {
      header.long.push({ refName: value, start: +start, end: +end })
    } else if (key === '#hole') {
      header.holes.push({ refName: value, start: +start, end: +end })
    } else {
      header.format ??= integerTag(line, 'format')
      header.chunk ??= integerTag(line, 'chunk')
      header.maxNode ??= integerTag(line, 'maxnode')
      header.revisited ??= integerTag(line, 'revisited')
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

// One stretch of a haplotype path: rows of consecutive piece index, decoded
export interface WalkRun {
  name: string
  fragStart: number
  hapOffset: number
  firstPiece: number
  // where each row after the first starts in `ids`, and the chunk each row
  // is filed under
  pieceStarts: number[]
  pieceChunks: number[]
  ids: Int32Array
  rev: Uint8Array
  // where each step lies on the reference, NaN off R or unknown, and whether
  // it runs against R's orientation at its pass
  refStarts: Float64Array
  against: Uint8Array
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
  const refStarts = new Float64Array(n).fill(Number.NaN)
  const against = new Uint8Array(n)
  const pieceStarts: number[] = []
  const pieceChunks: number[] = []
  let at = 0
  for (const p of pieces) {
    if (at > 0) {
      pieceStarts.push(at)
    }
    pieceChunks.push(p.chunkStart)
    ids.set(p.steps.ids, at)
    rev.set(p.steps.rev, at)
    if (p.steps.refStarts) {
      refStarts.set(p.steps.refStarts, at)
    }
    if (p.steps.against) {
      against.set(p.steps.against, at)
    }
    at += p.n
  }
  const first = pieces[0]!
  return {
    name: first.name,
    fragStart: first.fragStart,
    hapOffset: first.hapOffset,
    firstPiece: first.piece,
    pieceStarts,
    pieceChunks,
    ids,
    rev,
    refStarts,
    against,
    ...(first.abridged ? { abridged: first.abridged } : {}),
  }
}

/**
 * One walk's stretch through the cut, from its first step on the drawn
 * stretch of the reference to its last. `trimmedBefore` and `trimmedAfter`
 * count the bp of the walk the cut read but kept out, before and after it.
 */
export interface WalkFragment {
  name: string
  fragStart: number
  hapStart: number
  hapEnd: number
  ids: Int32Array
  rev: Uint8Array
  // where each step lies on the reference at its pass, NaN off it, and
  // whether it runs against the reference's orientation there
  refStarts: Float64Array
  against: Uint8Array
  trimmedBefore: number
  trimmedAfter: number
  // set on the fragment of an abridged row, from the piece's first step on R
  // to its last, which runs `abridged.bp` on the haplotype
  abridged?: Abridged
}

export interface WalkCutOptions {
  // keep the reference over the window and its context alone, and split each
  // walk where it steps onto the reference off them
  contained?: boolean
  // bp of each walk kept past either end of a contained span
  context?: number
  // the most bp of own sequence a departure gets its steps attached for
  attachBp?: number
  // a position with the reference's holes before it removed
  squeeze?: (bp: number) => number
}

// A record whose own steps the cut should attach and does not hold: the base
// of the reference beside its leave breakpoint names the chunk its leave
// piece is filed under, or the one before for a node straddling a boundary
export interface Unattached {
  departure: WalkDeparture
  base: number
}

/**
 * The cut a window draws from its walks: the reference inside the window plus
 * its context, the drawn stretch, and each run from its first step on that
 * stretch to its last. Nothing is followed past the drawn stretch and no
 * fragment runs on through nodes other walks kept, so a haplotype's fragment
 * does not depend on which others are cut.
 *
 * A contained cut splits a run wherever it steps onto the reference off the
 * drawn stretch. Each piece between those steps that steps on the stretch is
 * kept from its first such step to its last, so a detour that comes back is
 * kept whole, and `context` bp of the walk past either end. At AMY1
 * haplotypes go from a 4 kb window to the other copies and back: the cut
 * above keeps 9,928 nodes there, the contained cut 182.
 *
 * Every place a walk departs from the reference over the drawn stretch comes
 * back once: the departures inside one piece are computed from its steps,
 * and those between pieces, in different chunks or across W lines, come from
 * the file's records, which the reader read over the same stretch. A record
 * under `attachBp` of own sequence gets its own steps from the pieces the
 * cut holds, and the reader fetches the leave piece of one it does not.
 */
export function walkCut(
  runs: WalkRun[],
  nodes: Map<number, WalkNode>,
  refName: string,
  records: WalkDeparture[],
  {
    contained = false,
    context = 0,
    attachBp = ATTACH_BP,
    squeeze = (bp: number) => bp,
  }: WalkCutOptions = {},
) {
  // the window and its context, as the reference nodes on it span them
  let windowStart = Infinity
  let windowEnd = -Infinity
  for (const node of nodes.values()) {
    if (node.onWindow) {
      windowStart = Math.min(windowStart, node.start)
      windowEnd = Math.max(windowEnd, node.end)
    }
  }
  const length = (id: number) => {
    const node = nodes.get(id)
    return node === undefined ? 0 : node.end - node.start
  }
  const passesOf = (id: number) => nodes.get(id)?.passes ?? 0
  const onR = (run: WalkRun, i: number) => !Number.isNaN(run.refStarts[i])
  // bp of a run's steps before each, and after its last
  const prefixes = new Map<WalkRun, Float64Array>()
  const before = (run: WalkRun) => {
    let out = prefixes.get(run)
    if (!out) {
      out = new Float64Array(run.ids.length + 1)
      for (let i = 0; i < run.ids.length; i++) {
        out[i + 1] = out[i]! + length(run.ids[i]!)
      }
      prefixes.set(run, out)
    }
    return out
  }
  const stepOnWindow = (run: WalkRun, i: number) => {
    const start = run.refStarts[i]!
    const bp = before(run)
    return (
      !Number.isNaN(start) &&
      start < windowEnd &&
      start + bp[i + 1]! - bp[i]! > windowStart
    )
  }
  const spans: { run: WalkRun; first: number; last: number }[] = []
  // A contained run's spans: each stretch between steps onto the reference
  // off the window that has steps on it, from its first such step to its
  // last, and `context` bp of the walk past either
  const containedSpans = (run: WalkRun) => {
    const { ids } = run
    let from = 0
    let first = -1
    let last = -1
    for (let i = 0; i <= ids.length; i++) {
      const on = i < ids.length && stepOnWindow(run, i)
      if (i === ids.length || (onR(run, i) && !on)) {
        if (first >= 0) {
          for (let bp = 0; first > from && bp < context;) {
            bp += length(ids[--first]!)
          }
          for (let bp = 0; last < i - 1 && bp < context;) {
            bp += length(ids[++last]!)
          }
          spans.push({ run, first, last })
        }
        from = i + 1
        first = -1
      } else if (on) {
        if (first < 0) {
          first = i
        }
        last = i
      }
    }
  }
  const stretchSpan = (run: WalkRun) => {
    let first = -1
    let last = -1
    for (let i = 0; i < run.ids.length; i++) {
      if (stepOnWindow(run, i)) {
        if (first < 0) {
          first = i
        }
        last = i
      }
    }
    if (first >= 0) {
      spans.push({ run, first, last })
    }
  }
  // an abridged run is kept where its stretch of R overlaps the window
  const abridgedSpan = (run: WalkRun) => {
    const a = run.refStarts[0]!
    const b = run.refStarts[1]!
    if (Number.isNaN(a) || Number.isNaN(b)) {
      return
    }
    const ends = [a + length(run.ids[0]!), b + length(run.ids[1]!)]
    if (
      Math.max(ends[0]!, ends[1]!) > windowStart &&
      Math.min(a, b) < windowEnd
    ) {
      spans.push({ run, first: 0, last: 1 })
    }
  }
  for (const run of runs) {
    if (run.abridged) {
      abridgedSpan(run)
    } else if (contained) {
      containedSpans(run)
    } else {
      stretchSpan(run)
    }
  }
  const kept = new Set<number>()
  for (const [id, node] of nodes) {
    if (node.onWindow) {
      kept.add(id)
    }
  }
  const spansOf = new Map<WalkRun, { first: number; last: number }[]>()
  for (const span of spans) {
    for (let i = span.first; i <= span.last; i++) {
      kept.add(span.run.ids[i]!)
    }
    spansOf.set(span.run, [...(spansOf.get(span.run) ?? []), span])
  }
  const runEnd = (run: WalkRun) =>
    run.hapOffset + (run.abridged?.bp ?? before(run)[run.ids.length]!)
  const fragments: WalkFragment[] = []
  for (const { run, first, last } of spans) {
    const { ids, rev } = run
    if (run.abridged) {
      fragments.push({
        name: run.name,
        fragStart: run.fragStart,
        hapStart: run.hapOffset,
        hapEnd: runEnd(run),
        ids,
        rev,
        refStarts: run.refStarts,
        against: run.against,
        trimmedBefore: 0,
        trimmedAfter: 0,
        abridged: run.abridged,
      })
      continue
    }
    const bp = before(run)
    const hapStart = run.hapOffset + bp[first]!
    const hapEnd = run.hapOffset + bp[last + 1]!
    fragments.push({
      name: run.name,
      fragStart: run.fragStart,
      hapStart,
      hapEnd,
      ids: ids.subarray(first, last + 1),
      rev: rev.subarray(first, last + 1),
      refStarts: run.refStarts.subarray(first, last + 1),
      against: run.against.subarray(first, last + 1),
      trimmedBefore: hapStart - run.hapOffset,
      trimmedAfter: runEnd(run) - hapEnd,
    })
  }
  // A step's breakpoint on R: where the path's sequence stops following R,
  // the node's end for a step along R and its start for one against it, and
  // the other way round for a rejoin
  const breakpoint = (run: WalkRun, i: number, leave: boolean) => {
    const reversed = run.against[i] === 1
    const start = run.refStarts[i]!
    return {
      pos: reversed === leave ? start : start + length(run.ids[i]!),
      reversed,
    }
  }
  const overExtent = (d: WalkDeparture) => {
    const [lo, hi] = departureStretch(d)
    return lo < windowEnd && hi > windowStart
  }
  // where each walk's fragments start and end on its contig
  const fragmentEdges = new Map<string, Set<number>>()
  for (const f of fragments) {
    const edges = fragmentEdges.get(f.name) ?? new Set<number>()
    edges.add(f.hapStart)
    edges.add(f.hapEnd)
    fragmentEdges.set(f.name, edges)
  }
  const touchesFragment = (d: WalkDeparture) => {
    const edges = fragmentEdges.get(d.name)
    return (
      edges !== undefined &&
      (edges.has(d.hapStart) || edges.has(d.hapStart + d.ownBp))
    )
  }
  const departures: WalkDeparture[] = []
  for (const run of runs) {
    if (run.abridged || run.name === refName) {
      continue
    }
    const { ids, rev, refStarts, against } = run
    const bp = before(run)
    const bounds = pieceBounds(run)
    const spansHere = spansOf.get(run)
    const edges = fragmentEdges.get(run.name)
    const inSpan = (i: number, j: number) =>
      spansHere?.some(s => s.first <= i && j <= s.last) ?? false
    for (let p = 0; p + 1 < bounds.length; p++) {
      let last = -1
      for (let i = bounds[p]!; i < bounds[p + 1]!; i++) {
        const start = refStarts[i]!
        if (Number.isNaN(start)) {
          continue
        }
        if (last >= 0) {
          const ownBp = bp[i]! - bp[last + 1]!
          const leaveReversed = against[last] === 1
          const rejoinReversed = against[i] === 1
          const leavePos = leaveReversed
            ? refStarts[last]!
            : refStarts[last]! + bp[last + 1]! - bp[last]!
          const rejoinPos = rejoinReversed ? start + bp[i + 1]! - bp[i]! : start
          const adjacent =
            ownBp === 0 &&
            leaveReversed === rejoinReversed &&
            squeeze(leavePos) === squeeze(rejoinPos)
          if (!adjacent) {
            const leaveBase = leaveReversed ? leavePos : leavePos - 1
            const rejoinBase = rejoinReversed ? rejoinPos - 1 : rejoinPos
            const lo = Math.max(0, Math.min(leaveBase, rejoinBase))
            const hi = Math.max(lo, leaveBase, rejoinBase) + 1
            const hapStart = run.hapOffset + bp[last + 1]!
            if (
              (lo < windowEnd && hi > windowStart) ||
              inSpan(last, i) ||
              edges?.has(hapStart) === true ||
              edges?.has(hapStart + ownBp) === true
            ) {
              const d: WalkDeparture = {
                name: run.name,
                hapStart,
                ownBp,
                missingBp: 0,
                leavePos,
                leaveReversed,
                leavePasses: passesOf(ids[last]!),
                rejoinPos,
                rejoinReversed,
                rejoinPasses: passesOf(ids[i]!),
                fromSteps: true,
              }
              if (i > last + 1 && ownBp <= attachBp) {
                d.steps = {
                  ids: ids.subarray(last + 1, i),
                  rev: rev.subarray(last + 1, i),
                }
                for (let k = last + 1; k < i; k++) {
                  kept.add(ids[k]!)
                }
              }
              departures.push(d)
            }
          }
        }
        last = i
      }
    }
  }
  // A record's own steps: those after its leave step to the end of that
  // piece, then those before its rejoin step from the start of its piece. A
  // leave step is the last step on R of its piece and a rejoin step the
  // first of its, so the cut indexes those by name, breakpoint and contig
  // position
  const edgeKey = (name: string, pos: number, reversed: boolean, hap: number) =>
    `${name}\t${pos}\t${reversed ? '-' : '+'}\t${hap}`
  // a piece's last step on R and where the piece ends, and its first and
  // where the piece starts
  const leaves = new Map<string, [WalkRun, number, number]>()
  const rejoins = new Map<string, [WalkRun, number, number]>()
  for (const run of runs) {
    if (run.abridged) {
      continue
    }
    const bp = before(run)
    const bounds = pieceBounds(run)
    for (let p = 0; p + 1 < bounds.length; p++) {
      let first = -1
      let last = -1
      for (let i = bounds[p]!; i < bounds[p + 1]!; i++) {
        if (onR(run, i)) {
          if (first < 0) {
            first = i
          }
          last = i
        }
      }
      if (first < 0) {
        continue
      }
      const leave = breakpoint(run, last, true)
      leaves.set(
        edgeKey(
          run.name,
          leave.pos,
          leave.reversed,
          run.hapOffset + bp[last + 1]!,
        ),
        [run, last, bounds[p + 1]!],
      )
      const rejoin = breakpoint(run, first, false)
      rejoins.set(
        edgeKey(
          run.name,
          rejoin.pos,
          rejoin.reversed,
          run.hapOffset + bp[first]!,
        ),
        [run, first, bounds[p]!],
      )
    }
  }
  const unattached: Unattached[] = []
  for (const d of records) {
    if (d.ownBp > d.missingBp && d.ownBp <= attachBp) {
      const ids: number[] = []
      const rev: number[] = []
      const leave = leaves.get(
        edgeKey(d.name, d.leavePos, d.leaveReversed, d.hapStart),
      )
      if (leave) {
        const [run, i, end] = leave
        for (let k = i + 1; k < end; k++) {
          ids.push(run.ids[k]!)
          rev.push(run.rev[k]!)
        }
      }
      const rejoin = rejoins.get(
        edgeKey(d.name, d.rejoinPos, d.rejoinReversed, d.hapStart + d.ownBp),
      )
      if (rejoin) {
        const [run, i, start] = rejoin
        for (let k = start; k < i; k++) {
          ids.push(run.ids[k]!)
          rev.push(run.rev[k]!)
        }
      }
      if (ids.length > 0) {
        d.steps = { ids: Int32Array.from(ids), rev: Uint8Array.from(rev) }
        for (const id of ids) {
          kept.add(id)
        }
      }
      if (!leave) {
        unattached.push({
          departure: d,
          base: d.leaveReversed ? d.leavePos : Math.max(0, d.leavePos - 1),
        })
      }
    }
    if (overExtent(d) || touchesFragment(d)) {
      departures.push(d)
    }
  }
  departures.sort(
    (a, b) =>
      (a.name < b.name ? -1 : a.name > b.name ? 1 : 0) ||
      a.hapStart - b.hapStart ||
      a.leavePos - b.leavePos,
  )
  return { kept, fragments, departures, unattached }
}

// Where each of a run's pieces starts, and where the last ends
function pieceBounds(run: WalkRun) {
  const bounds = [0]
  for (const start of run.pieceStarts) {
    bounds.push(start)
  }
  bounds.push(run.ids.length)
  return bounds
}

/**
 * Adds to `into` the nodes of decoded rows that walkCut keeps whatever else
 * the cut reads: those of `chunk.nodes` on the window, and those a row steps
 * on between two of its steps on the window, `start`-`end` with its context,
 * with no step onto the reference off the window between them in a contained
 * cut. A node placed later at a reference position off the window can split a
 * contained walk, so the count holds for the placements at the time.
 */
export function surelyKept(
  chunk: { nodes: Map<number, unknown>; rows: WalkRow[] },
  nodes: Map<number, WalkNode>,
  { contained, start, end }: { contained: boolean; start: number; end: number },
  into: Set<number>,
) {
  for (const id of chunk.nodes.keys()) {
    if (nodes.get(id)?.onWindow) {
      into.add(id)
    }
  }
  for (const { steps } of chunk.rows) {
    const { ids, refStarts } = steps
    let last = -1
    for (let i = 0; i < ids.length; i++) {
      const node = nodes.get(ids[i]!)
      const at = refStarts?.[i]
      const onR = at !== undefined && !Number.isNaN(at)
      if (
        node?.onWindow &&
        onR &&
        at < end &&
        at + node.end - node.start > start
      ) {
        for (let j = last + 1; last >= 0 && j < i; j++) {
          into.add(ids[j]!)
        }
        last = i
      } else if (contained && onR) {
        last = -1
      }
    }
  }
}

// Whether placement `node` beats `held`: a pass of the reference over the
// window beats any other, then a reference position beats one off it, and
// among equals the path sorting first at its earliest start wins, so the order
// a cut reads chunks in does not matter
const outranks = (node: WalkNode, held: WalkNode) =>
  node.onWindow !== held.onWindow
    ? node.onWindow
    : node.rank === 0 &&
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
    const { window } = this
    const placed = {
      ...node,
      onWindow:
        node.rank === 0 &&
        node.refName === window.refName &&
        node.end > window.start - window.context &&
        node.start < window.end + window.context,
    }
    const held = this.nodes.get(id)
    if (!held || outranks(placed, held)) {
      this.nodes.set(id, placed)
    }
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
   * name, and the departures with their own steps. graphTablesGFA writes
   * them as the GFA they stand for, the same bytes for the same rows.
   */
  tables(
    kept: Set<number>,
    fragments: WalkFragment[],
    departures: WalkDeparture[] = [],
  ): GraphTables {
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
    const passes = new Int32Array(held.length)
    held.forEach(({ id, refName, start, end, rank, passes: n }, i) => {
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
      passes[i] = n
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
      nodes: { names, lengths, refs, starts, ranks, refNames, passes },
      links: { from, to, strands },
      walks: {
        names: walks.map(f => f.name),
        starts: Float64Array.from(walks, f => f.hapStart),
        ends: Float64Array.from(walks, f => f.hapEnd),
        offsets,
        steps,
        reversed,
        trimmedBefore: Float64Array.from(walks, f => f.trimmedBefore),
        trimmedAfter: Float64Array.from(walks, f => f.trimmedAfter),
        abridgedSteps: Float64Array.from(walks, f => f.abridged?.steps ?? NaN),
        abridgedBp: Float64Array.from(walks, f => f.abridged?.bp ?? NaN),
      },
      departures: departureTable(departures, nodeIndex),
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

function departureTable(
  departures: WalkDeparture[],
  nodeIndex: (id: number) => number,
): DepartureTable {
  const n = departures.length
  const offsets = new Int32Array(n + 1)
  departures.forEach((d, k) => {
    offsets[k + 1] = offsets[k]! + (d.steps?.ids.length ?? 0)
  })
  const steps = new Int32Array(offsets[n]!)
  const reversed = new Uint8Array(steps.length)
  departures.forEach((d, k) => {
    if (d.steps) {
      const at = offsets[k]!
      for (let j = 0; j < d.steps.ids.length; j++) {
        steps[at + j] = nodeIndex(d.steps.ids[j]!)
      }
      reversed.set(d.steps.rev, at)
    }
  })
  return {
    names: departures.map(d => d.name),
    hapStarts: Float64Array.from(departures, d => d.hapStart),
    ownBp: Float64Array.from(departures, d => d.ownBp),
    missingBp: Float64Array.from(departures, d => d.missingBp),
    leavePos: Float64Array.from(departures, d => d.leavePos),
    leaveReversed: Uint8Array.from(departures, d => (d.leaveReversed ? 1 : 0)),
    leavePasses: Int32Array.from(departures, d => d.leavePasses),
    rejoinPos: Float64Array.from(departures, d => d.rejoinPos),
    rejoinReversed: Uint8Array.from(departures, d =>
      d.rejoinReversed ? 1 : 0,
    ),
    rejoinPasses: Int32Array.from(departures, d => d.rejoinPasses),
    fromSteps: Uint8Array.from(departures, d => (d.fromSteps ? 1 : 0)),
    offsets,
    steps,
    reversed,
  }
}

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

/**
 * The zoom-in notice for a cut whose walk rows hold more steps than
 * `readPerKept` times `stepBudget`, or undefined when they fit. A window reads
 * whole chunks, those before it included, so the span that fits is counted in
 * chunks from the window's start; when the chunks a window there always reads
 * are over on their own, no zoom fits and the notice asks for fewer
 * haplotypes. A cut for every haplotype is also offered fewer where the
 * reference and one haplotype, at the walks' average steps, would fit the
 * whole window.
 */
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
 * The zoom-in notice for a cut that keeps more steps than `budget`, the
 * fragments' and the departures' attached steps together, or undefined when
 * it fits. The span that fits is the window scaled by how far over the cut
 * is.
 */
export function keptStepsError(
  fragments: WalkFragment[],
  departures: WalkDeparture[],
  budget: number,
  window: { start: number; end: number },
  filtered: boolean,
) {
  const total =
    fragments.reduce((sum, f) => sum + f.ids.length, 0) +
    departures.reduce((sum, d) => sum + (d.steps?.ids.length ?? 0), 0)
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
 * The zoom-in notice for a cut whose walk read would fetch more than
 * `budget` compressed bytes, or undefined when it fits. `bytesTo(end)` is the
 * index's estimate for a read from the window's first query base to `end`;
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
        { bytes, zoomHelps: true },
      )
    : new RegionTooLargeError(
        `Too much graph here to fetch, even zoomed in (${megabytes(least ?? bytes)} against ${names.bytes} ${megabytes(budget)}). A satellite or rDNA array usually causes this: the graph folds the repeat into a few nodes that each haplotype loops through many times`,
        budget,
        windowBp,
        1,
        { bytes },
      )
}

export const megabytes = (bytes: number) => `${(bytes / 1e6).toFixed(1)} MB`
