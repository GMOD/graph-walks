import { EARLIER_FORMAT_6, WalkRowError } from './walkChunk.ts'

// A departure: where a path's consecutive steps on the reference sequence S
// are not adjacent there, from a walk file's `departure` records or computed
// from a piece's steps (gfa-to-tabix's docs/walks-format.md, ADR 0009). Each
// end is the step's breakpoint on S, where the path's sequence stops
// following S, whether the step runs against S's orientation at its pass, and
// how many times the reference passes the step's node.

export interface WalkDeparture {
  name: string
  // where the path's own sequence between the ends starts on its contig, right
  // after the leave step, and how far it runs: 0 for a direct jump
  hapStart: number
  ownBp: number
  // the part of ownBp between the contig's W lines, which the graph lacks
  missingBp: number
  leavePos: number
  leaveReversed: boolean
  // how many times the reference's paths pass the end's node: above 1, the
  // pass is the writer's placement choice, and so is the departure
  leavePasses: number
  rejoinPos: number
  rejoinReversed: boolean
  rejoinPasses: number
  // computed from a piece's steps by the cut, rather than read from a record
  fromSteps: boolean
  // the own steps the cut holds, in the path's order; fewer bp than ownBp
  // less missingBp where some lie on other reference sequences, or in
  // pieces the cut did not read
  steps?: { ids: Int32Array; rev: Uint8Array }
}

// A stretch of a reference sequence no line of the reference covers, from a
// `#hole` header line
export interface Hole {
  refName: string
  start: number
  end: number
}

// The most bp of own sequence a departure gets its steps attached for
export const ATTACH_BP = 65_536

// A position on `refName` with the holes before it removed, so that steps
// either side of a hole compare as adjacent
export function squeezer(holes: Hole[], refName: string) {
  const spans = holes
    .filter(h => h.refName === refName)
    .sort((a, b) => a.start - b.start)
  return spans.length === 0
    ? (bp: number) => bp
    : (bp: number) => {
        let out = bp
        for (const { start, end } of spans) {
          if (end > bp) {
            break
          }
          out -= end - start
        }
        return out
      }
}

type Ends = Pick<
  WalkDeparture,
  'ownBp' | 'leavePos' | 'leaveReversed' | 'rejoinPos' | 'rejoinReversed'
>

// Whether two consecutive steps on S follow S at their passes, and so depart
// nowhere: no sequence of the path's own between them, the same strand and
// one's breakpoint the other's, holes removed
export function adjacentSteps(d: Ends, squeeze: (bp: number) => number) {
  return (
    d.ownBp === 0 &&
    d.leaveReversed === d.rejoinReversed &&
    squeeze(d.leavePos) === squeeze(d.rejoinPos)
  )
}

/**
 * What a departure does on S, derived from its ends: `switch`, the strands
 * differ; `back`, the path jumps back over reference it walks again;
 * `adjacent`, sequence of its own between adjacent steps; else `forward`, a
 * jump forward over reference it leaves out. Positions compare with the
 * holes before them removed.
 */
export type DepartureKind = 'switch' | 'back' | 'adjacent' | 'forward'

export function departureKind(
  d: Ends,
  squeeze: (bp: number) => number = bp => bp,
): DepartureKind {
  if (d.leaveReversed !== d.rejoinReversed) {
    return 'switch'
  }
  const leave = squeeze(d.leavePos)
  const rejoin = squeeze(d.rejoinPos)
  return leave === rejoin
    ? 'adjacent'
    : (d.leaveReversed ? rejoin > leave : rejoin < leave)
      ? 'back'
      : 'forward'
}

// The bases of S a departure covers, from the base beside one breakpoint to
// the base beside the other, [lo, hi), as its record is filed
export function departureStretch(d: Ends): [number, number] {
  const leaveBase = d.leaveReversed ? d.leavePos : d.leavePos - 1
  const rejoinBase = d.rejoinReversed ? d.rejoinPos - 1 : d.rejoinPos
  const lo = Math.max(0, Math.min(leaveBase, rejoinBase))
  return [lo, Math.max(lo, leaveBase, rejoinBase) + 1]
}

/**
 * Whether a departure's place on the reference is the writer's choice: an
 * end lies on a node the reference passes more than once, so another pass
 * would put it elsewhere or make it none, as where a path follows one copy
 * of a collapsed repeat and the placement hops between copies
 */
export function placementChoice(
  d: Pick<WalkDeparture, 'leavePasses' | 'rejoinPasses'>,
) {
  return d.leavePasses > 1 || d.rejoinPasses > 1
}

const count = (field: string | undefined, row: string[], least = 0) => {
  const n = Number(field)
  if (
    field === undefined ||
    field === '' ||
    !Number.isInteger(n) ||
    n < least
  ) {
    throw new WalkRowError(`not a walk format 6 row: ${row.join('\t')}`)
  }
  return n
}

const strand = (field: string | undefined, row: string[]) => {
  if (field !== '+' && field !== '-') {
    throw new WalkRowError(`not a walk format 6 row: ${row.join('\t')}`)
  }
  return field === '-'
}

// A record row's columns, `departure path hapStart ownBp missingBp leavePos
// leaveStrand leavePasses rejoinPos rejoinStrand rejoinPasses`, from its kind
// at `from`; a row that is not one throws, saying to rebuild a file whose
// rows lack the pass counts
export function parseDepartureRow(cols: string[], from = 3): WalkDeparture {
  const n = (k: number) => count(cols[from + k], cols)
  const passes = (k: number) => count(cols[from + k], cols, 1)
  if (cols.length - from === 9) {
    throw new WalkRowError(
      `a departure record with no pass counts, ${EARLIER_FORMAT_6}: ${cols.join('\t')}`,
    )
  }
  return {
    name: cols[from + 1]!,
    hapStart: n(2),
    ownBp: n(3),
    missingBp: n(4),
    leavePos: n(5),
    leaveReversed: strand(cols[from + 6], cols),
    leavePasses: passes(7),
    rejoinPos: n(8),
    rejoinReversed: strand(cols[from + 9], cols),
    rejoinPasses: passes(10),
    fromSteps: false,
  }
}

// A record row's `ownTiles`, after `rejoinPasses`: the nominal starts of
// the tiles holding its own steps, `.` for none
export function departureOwnTiles(cols: string[], from = 3) {
  const field = cols[from + 11]
  const starts = field === '.' || field === undefined ? [] : field.split(',')
  return starts.map(start => {
    const n = Number(start)
    if (start === '' || !Number.isSafeInteger(n) || n < 0) {
      throw new WalkRowError(
        `a departure record's ownTiles is not a list of tile starts: ${cols.join('\t')}`,
      )
    }
    return n
  })
}

// Every record's rows repeat the same columns after the block, so this keys
// a row's record
export function departureRowKey(cols: string[], from = 3) {
  return cols.slice(from).join('\t')
}
