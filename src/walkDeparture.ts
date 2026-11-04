// A departure: where a path's consecutive steps on the reference sequence S
// are not adjacent there, from a walk file's `departure` records or computed
// from a piece's steps (gfa-to-tabix's docs/walks-format.md, ADR 0009). Each
// end is the step's breakpoint on S, where the path's sequence stops
// following S, and whether the step runs against S's orientation at its pass.

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
  rejoinPos: number
  rejoinReversed: boolean
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

// A record row's columns, `departure path hapStart ownBp missingBp leavePos
// leaveStrand rejoinPos rejoinStrand`, from its kind at `from`
export function parseDepartureRow(cols: string[], from = 3): WalkDeparture {
  return {
    name: cols[from + 1]!,
    hapStart: +cols[from + 2]!,
    ownBp: +cols[from + 3]!,
    missingBp: +cols[from + 4]!,
    leavePos: +cols[from + 5]!,
    leaveReversed: cols[from + 6] === '-',
    rejoinPos: +cols[from + 7]!,
    rejoinReversed: cols[from + 8] === '-',
    fromSteps: false,
  }
}

// Every record's rows repeat the same columns after the block, so this keys
// a row's record
export function departureRowKey(cols: string[], from = 3) {
  return cols.slice(from).join('\t')
}
