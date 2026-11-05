import {
  RegionTooLargeError,
  bpString,
  panSNHaplotype,
  panSNSample,
} from './util.ts'

import type { WalkDeparture } from './walkDeparture.ts'
import type { WalkReader, WalkReaderHeader } from './walkReader.ts'
import type { WalkFragment, WalkNode } from './walkRows.ts'

/**
 * Where a haplotype leaves the reference path R and rejoins it: R's bases
 * `refStart`–`refEnd` against `altLength` of its own. `invertedBp` counts the
 * steps along R against its orientation. A step back along R walks a tandem
 * copy again: the departure is then a point at the junction, `repeatedBp`
 * counts the bp walked again, and `copied` spans them on R. An `abridged`
 * departure stands for the pieces a writer abridged, so only its net change
 * is exact.
 */
export interface Departure {
  refStart: number
  refEnd: number
  altLength: number
  invertedBp: number
  repeatedBp: number
  copied?: [number, number]
  abridged?: true
}

/**
 * A `sample#haplotype`'s departures, by refStart. `covered` spans the R
 * exactly one of its paths reaches, each span within one path.
 */
export interface HaplotypeDepartures {
  haplotype: string
  covered: [number, number][]
  departures: Departure[]
}

// A cut's fragments, its departures and its nodes, as WalkReader.cut returns
// them
export interface DepartureSource {
  fragments: WalkFragment[]
  departures: WalkDeparture[]
  nodes: Map<number, WalkNode>
}

// What a path does along its contig: follows R over [lo, hi), along or
// against it; runs `bp` of its own; or runs `bp` the writer abridged between
// R's bases lo and hi, along R or against it as its ends are
type Event =
  | { kind: 'on'; lo: number; hi: number; reversed: boolean }
  | { kind: 'own'; bp: number }
  | { kind: 'abridged'; lo: number; hi: number; bp: number; reversed: boolean }

// A point on a path's contig: where it is on R until (a fragment's last step
// or a departure's leave, `order` 0), from (a fragment's first step or a
// departure's rejoin, 1), or where an abridged stretch starts (2). `own` is
// the bp of the path's own that follow a leave
interface Anchor {
  hap: number
  order: 0 | 1 | 2
  pos: number
  reversed: boolean
  exit?: true
  own?: number
  abridged?: { lo: number; hi: number; bp: number; reversed: boolean }
}

// A fragment or departure as the bp of its contig it holds, the R it reaches
// and the bp of the path it stands for
type Item = {
  hap: [number, number]
  lo: number
  hi: number
  bp: number
} & ({ fragment: WalkFragment } | { departure: WalkDeparture })

const onR = (start: number) => !Number.isNaN(start)
// the lower and the higher of two positions, those off R left out
const least = (a: number, b: number) =>
  onR(a) ? (onR(b) ? Math.min(a, b) : a) : b
const most = (a: number, b: number) =>
  onR(a) ? (onR(b) ? Math.max(a, b) : a) : b

/**
 * Each haplotype's departures from R, `refName`, read off a cut: its
 * departures, every one over the drawn stretch, give where each path leaves
 * R and rejoins it and how much of its own it runs between, and its
 * fragments' first and last steps give where it is on R otherwise. A
 * fragment with the departures at its ends is a chain, as is a departure
 * touching no fragment, such as one jumping over the window; each reads in
 * its own orientation, mostly reversed reading reversed, and a stretch of R
 * walked against it reads as inverted, a jump back as a tandem copy, and
 * sequence between a contig's W lines as the path's own. R that two chains
 * of a haplotype reach, such as a contig's second visit to the window, has
 * no call. An abridged fragment's bp, less the own bp of the departures
 * leaving or rejoining inside it, are a net change between the R its handles
 * reach, along R or against it as both handles are. The walks on R's own
 * haplotype are left out.
 */
export function walkDepartures(
  { fragments, departures, nodes }: DepartureSource,
  refName: string,
) {
  const length = (id: number) => {
    const node = nodes.get(id)
    return node ? node.end - node.start : 0
  }
  const reference = panSNHaplotype(refName) ?? refName
  const byPath = new Map<string, Item[]>()
  const add = (name: string, item: Item) => {
    const items = byPath.get(name)
    if (items) {
      items.push(item)
    } else {
      byPath.set(name, [item])
    }
  }
  for (const f of fragments) {
    add(f.name, fragmentItem(f, length))
  }
  for (const d of departures) {
    add(d.name, {
      hap: [d.hapStart, d.hapStart + d.ownBp],
      lo: Math.min(d.leavePos, d.rejoinPos),
      hi: Math.max(d.leavePos, d.rejoinPos),
      bp: d.ownBp,
      departure: d,
    })
  }
  const byHaplotype = new Map<string, HaplotypeDepartures>()
  for (const [name, items] of byPath) {
    const haplotype = panSNHaplotype(name) ?? panSNSample(name)
    if (haplotype === reference || name === refName) {
      continue
    }
    let h = byHaplotype.get(haplotype)
    if (!h) {
      h = { haplotype, covered: [], departures: [] }
      byHaplotype.set(haplotype, h)
    }
    for (const chain of pathChains(items, length)) {
      readPath(chain, h)
    }
  }
  const out = [...byHaplotype.values()]
  for (const h of out) {
    h.departures.sort((a, b) => a.refStart - b.refStart || a.refEnd - b.refEnd)
    h.covered = coveredOnce(h.covered)
  }
  return out
}

function fragmentItem(f: WalkFragment, length: (id: number) => number): Item {
  const n = f.ids.length
  const sa = f.refStarts[0]!
  const sb = f.refStarts[n - 1]!
  const lo = least(sa, sb)
  const hi = most(sa + length(f.ids[0]!), sb + length(f.ids[n - 1]!))
  return f.abridged
    ? { hap: [f.hapStart, f.hapEnd], lo, hi, bp: f.abridged.bp, fragment: f }
    : {
        hap: [f.hapStart - f.trimmedBefore, f.hapEnd + f.trimmedAfter],
        lo,
        hi,
        bp: f.hapEnd + f.trimmedAfter - f.hapStart + f.trimmedBefore,
        fragment: f,
      }
}

/**
 * A path's chains, each its events in the path's order: a fragment with the
 * departures leaving its last step, rejoining its first or lying inside it,
 * fragments that abut, and a departure touching no fragment on its own. What
 * the read did not see lies between chains, so each reads on its own, and R
 * that two chains reach has no call. An abridged fragment is a net over one
 * piece a writer moved out of a chunk, from its row's hapOffset to hapOffset
 * + bp on the contig. A departure keeps all its own bp, so a piece it leaves
 * inside ends at its leave and one it rejoins inside starts at its rejoin:
 * those steps of the piece are the departure's own.
 */
function pathChains(items: Item[], length: (id: number) => number) {
  const departures = items.filter(i => 'departure' in i)
  const pieces = items.flatMap(item =>
    'fragment' in item && item.fragment.abridged
      ? [ownTrimmed(item, departures)]
      : [],
  )
  const fragments = items.filter(i => 'fragment' in i && !i.fragment.abridged)
  // where each item starts and ends on the contig
  const edges = (item: Item): [number, number] =>
    'fragment' in item && !item.fragment.abridged
      ? [item.fragment.hapStart, item.fragment.hapEnd]
      : item.hap
  // items chain where one's contig stretch touches another's: a departure's
  // leave at a fragment's last step or inside it, its rejoin at the first,
  // or two fragments abutting. In start order, an item touches the chain so
  // far when it starts no later than the chain's furthest end
  const all = [...pieces, ...fragments, ...departures]
    .map(item => ({ item, edges: edges(item) }))
    .sort((a, b) => a.edges[0] - b.edges[0])
  const chains: Item[][] = []
  let chainEnd = -Infinity
  for (const {
    item,
    edges: [start, end],
  } of all) {
    const chain = chains.at(-1)
    if (chain && start <= chainEnd) {
      chain.push(item)
    } else {
      chains.push([item])
    }
    chainEnd = Math.max(chainEnd, end)
  }
  return chains.map(chain => chainEvents(chain, length))
}

// An abridged piece without the own steps of the departures that leave or
// rejoin inside it: its steps on R lie outside every departure's own
// sequence, so those are its leading and trailing steps
function ownTrimmed(piece: Item, departures: Item[]): Item {
  let [start, end] = piece.hap
  for (const { hap } of departures) {
    if (hap[0] > start && hap[0] < end) {
      end = hap[0]
    }
    if (hap[1] > start && hap[1] < end) {
      start = hap[1]
    }
  }
  return {
    ...piece,
    hap: [start, end],
    bp: piece.bp - (start - piece.hap[0]) - (piece.hap[1] - end),
  }
}

function chainEvents(items: Item[], length: (id: number) => number) {
  const anchors: Anchor[] = []
  for (const item of items) {
    if ('fragment' in item) {
      const f = item.fragment
      if (f.abridged) {
        if (onR(item.lo) && onR(item.hi)) {
          const reversed = f.against.every(r => r === 1)
          anchors.push({
            hap: item.hap[0],
            order: 2,
            pos: item.lo,
            reversed,
            abridged: { lo: item.lo, hi: item.hi, bp: item.bp, reversed },
          })
        }
        continue
      }
      const n = f.ids.length
      const sa = f.refStarts[0]!
      const sb = f.refStarts[n - 1]!
      const ra = f.against[0] === 1
      const rb = f.against[n - 1] === 1
      if (onR(sa)) {
        const pos = ra ? sa + length(f.ids[0]!) : sa
        anchors.push({ hap: f.hapStart, order: 1, pos, reversed: ra })
      }
      if (onR(sb)) {
        const pos = rb ? sb : sb + length(f.ids[n - 1]!)
        anchors.push({ hap: f.hapEnd, order: 0, pos, reversed: rb, exit: true })
      }
    } else {
      const d = item.departure
      anchors.push(
        {
          hap: d.hapStart,
          order: 0,
          pos: d.leavePos,
          reversed: d.leaveReversed,
          own: d.ownBp,
        },
        {
          hap: d.hapStart + d.ownBp,
          order: 1,
          pos: d.rejoinPos,
          reversed: d.rejoinReversed,
        },
      )
    }
  }
  anchors.sort((a, b) => a.hap - b.hap || a.order - b.order)
  const events: Event[] = []
  const on = (from: number, to: number, reversed: boolean): Event => ({
    kind: 'on',
    lo: Math.min(from, to),
    hi: Math.max(from, to),
    reversed,
  })
  // the open stretch on R: where it began, and the furthest fragment end seen
  // along it
  let cur: { pos: number; reversed: boolean } | undefined
  let end: number | undefined
  const close = () => {
    if (cur) {
      events.push(on(cur.pos, end ?? cur.pos, cur.reversed))
    }
    cur = undefined
    end = undefined
  }
  for (const a of anchors) {
    if (cur && !a.abridged && cur.reversed !== a.reversed) {
      close()
    }
    if (a.abridged) {
      close()
      events.push({ kind: 'abridged', ...a.abridged })
    } else if (a.order === 1) {
      cur ??= { pos: a.pos, reversed: a.reversed }
    } else if (a.exit) {
      if (cur) {
        end = a.pos
      }
    } else {
      cur ??= { pos: a.pos, reversed: a.reversed }
      end = a.pos
      close()
      events.push({ kind: 'own', bp: a.own ?? 0 })
    }
  }
  close()
  return events
}

/**
 * A path's departures in the forward view: its stretches along R in order,
 * reversed where it walks R mostly against its orientation, each stretch
 * against it or of its own between two along it making one departure.
 */
function readPath(events: Event[], h: HaplotypeDepartures) {
  let bias = 0
  for (const e of events) {
    if (e.kind !== 'own') {
      bias += (e.reversed ? 1 : -1) * (e.hi - e.lo || 1)
    }
  }
  const flip = bias > 0
  let first = NaN
  let prevEnd = NaN
  let reach = -Infinity
  let alt = 0
  let inverted = 0
  // set while an abridged stretch's bp are in `alt`: the far end of the R
  // its handles reach where it runs along R, and true where it runs against
  let abridgedEnd = NaN
  let abridgedAgainst = false
  for (const e of flip ? events.toReversed() : events) {
    if (e.kind === 'own') {
      alt += e.bp
    } else if (e.kind === 'abridged' && e.reversed !== flip) {
      alt += e.bp
      inverted += e.hi - e.lo
      abridgedAgainst = true
    } else if (e.kind === 'abridged') {
      if (Number.isNaN(prevEnd)) {
        first = e.lo
        prevEnd = e.lo
      }
      alt += e.bp
      abridgedEnd = Number.isNaN(abridgedEnd)
        ? e.hi
        : Math.max(abridgedEnd, e.hi)
    } else if (e.reversed !== flip) {
      alt += e.hi - e.lo
      inverted += e.hi - e.lo
    } else {
      const start = e.lo
      if (Number.isNaN(prevEnd)) {
        first = start
      } else if (!Number.isNaN(abridgedEnd)) {
        h.departures.push({
          refStart: prevEnd,
          refEnd: Math.max(start, prevEnd),
          altLength: alt + Math.max(0, prevEnd - start),
          invertedBp: inverted,
          repeatedBp: 0,
          abridged: true,
        })
        reach = Math.max(reach, abridgedEnd)
        abridgedEnd = NaN
      } else if (start !== prevEnd || alt > 0) {
        const repeated = start < prevEnd ? prevEnd - start : 0
        h.departures.push({
          refStart: prevEnd,
          refEnd: repeated > 0 ? prevEnd : start,
          altLength: alt + repeated,
          invertedBp: inverted,
          repeatedBp: repeated,
          ...(repeated > 0 ? { copied: [start, prevEnd] } : {}),
          ...(abridgedAgainst ? { abridged: true as const } : {}),
        })
      }
      prevEnd = e.hi
      reach = Math.max(reach, e.hi)
      alt = 0
      inverted = 0
      abridgedAgainst = false
    }
  }
  if (!Number.isNaN(abridgedEnd)) {
    h.departures.push({
      refStart: prevEnd,
      refEnd: Math.max(abridgedEnd, prevEnd),
      altLength: alt,
      invertedBp: inverted,
      repeatedBp: 0,
      abridged: true,
    })
    reach = Math.max(reach, abridgedEnd)
  }
  if (!Number.isNaN(first)) {
    h.covered.push([first, reach])
  }
}

/**
 * The stretches of R that exactly one of a haplotype's paths covers, each
 * within one path. Where two paths overlap, such as a second contig over the
 * same locus, the haplotype holds two copies of that stretch and no call;
 * where one path ends and the next starts, no call crosses the junction.
 */
function coveredOnce(spans: [number, number][]) {
  const points = [...new Set(spans.flat())].sort((a, b) => a - b)
  const out: [number, number][] = []
  let lastPath = -1
  for (let k = 0; k + 1 < points.length; k++) {
    const a = points[k]!
    const b = points[k + 1]!
    const paths = spans.flatMap(([s, e], i) => (s <= a && e >= b ? [i] : []))
    const last = out.at(-1)
    if (paths.length !== 1) {
      lastPath = -1
    } else if (last && paths[0] === lastPath && last[1] === a) {
      last[1] = b
    } else {
      out.push([a, b])
      lastPath = paths[0]!
    }
  }
  return out
}

// The most steps a forms cut keeps, expanded abridged walks included; an
// expanded cut at the MHC keeps 9.4 M
export const FORMS_STEP_BUDGET = 16_000_000

export interface FormsCutOptions {
  // read abridged walks in full, as far as the step budget allows; true by
  // default
  expand?: boolean
  // false reads every abridged walk in full, or throws ExpansionTooLargeError
  // before the second cut where that takes the cut past stepBudget; true by
  // default, which reads them fewest steps first while the budget lasts
  partial?: boolean
  // the most steps the cut keeps with the abridged walks it reads in full,
  // FORMS_STEP_BUDGET by default; each cut keeps at least FORMS_STEP_BUDGET
  stepBudget?: number
  // the most compressed bytes either cut fetches from the region's chunks,
  // past which the cut throws RegionTooLargeError; WALK_BUDGETS' by default.
  // Only stepBudget bounds the detail records expansion reads
  byteBudget?: number
  // the most compressed bytes either cut's departure read fetches, past
  // which the cut throws DepartureReadError; WALK_BUDGETS' by default
  departureByteBudget?: number
  signal?: AbortSignal
}

/**
 * A forms cut with `partial: false` whose abridged walks, read in full, would
 * keep `steps` past `limit`, the step budget. `fitsBp` scales the region by
 * how far over it is.
 */
export class ExpansionTooLargeError extends RegionTooLargeError {
  readonly steps: number

  constructor(steps: number, limit: number, windowBp: number) {
    const fitsBp = Math.max(1, Math.floor((windowBp * limit) / steps))
    super(
      `Zoom in to about ${bpString(fitsBp)} to read every walk in full`,
      limit,
      windowBp,
      fitsBp,
    )
    this.steps = steps
  }
}

/**
 * Every haplotype's departures over a region, for structuralForms, from a
 * cut without links: its departures hold every one over the region, those
 * leaving it and rejoining far away or jumping over it whole included.
 * Abridged walks are read in full, those with the fewest steps first, while
 * the cut stays under `stepBudget`; any left keep their net change only.
 * `partial: false` reads every abridged walk or throws
 * ExpansionTooLargeError. A haplotype the header names with no walk in the
 * cut is listed with nothing covered. A region the index estimates past
 * `byteBudget` throws the reader's RegionTooLargeError before any read.
 */
export async function formsCut(
  reader: Pick<WalkReader, 'header' | 'sequenceFor' | 'cut'>,
  region: { refName: string; start: number; end: number },
  opts: FormsCutOptions = {},
) {
  const {
    expand = true,
    partial = true,
    stepBudget = FORMS_STEP_BUDGET,
    byteBudget,
    departureByteBudget,
    signal,
  } = opts
  const header: WalkReaderHeader = await reader.header({ signal })
  const sample = header.references[0] ?? ''
  const refName =
    (await reader.sequenceFor(sample, region.refName, { signal })) ??
    region.refName
  const cutOpts = {
    links: false,
    stepBudget: Math.max(stepBudget, FORMS_STEP_BUDGET),
    byteBudget,
    departureByteBudget,
    signal,
  }
  let cut = await reader.cut(refName, region, cutOpts)
  if (expand) {
    const expanded = expandedWalks(cut.fragments, stepBudget, partial, region)
    if (expanded.size > 0) {
      cut = await reader.cut(refName, region, {
        ...cutOpts,
        expand: name => expanded.has(name),
      })
    }
  }
  const haplotypes = walkDepartures(
    {
      fragments: cut.fragments,
      departures: cut.departures,
      nodes: cut.graph.nodes,
    },
    refName,
  )
  const seen = new Set(haplotypes.map(h => h.haplotype))
  for (const haplotype of header.haplotypes) {
    if (!seen.has(haplotype)) {
      haplotypes.push({ haplotype, covered: [], departures: [] })
    }
  }
  return haplotypes
}

// The abridged walks a cut can read in full within `stepBudget`, fewest
// steps first, or all of them where `partial` is false
function expandedWalks(
  fragments: WalkFragment[],
  stepBudget: number,
  partial: boolean,
  region: { start: number; end: number },
) {
  let kept = 0
  const abridged = new Map<string, number>()
  for (const f of fragments) {
    if (f.abridged) {
      abridged.set(f.name, (abridged.get(f.name) ?? 0) + f.abridged.steps)
    } else {
      kept += f.ids.length
    }
  }
  if (!partial) {
    const steps = [...abridged.values()].reduce((sum, n) => sum + n, kept)
    if (abridged.size > 0 && steps > stepBudget) {
      throw new ExpansionTooLargeError(
        steps,
        stepBudget,
        region.end - region.start,
      )
    }
    return new Set(abridged.keys())
  }
  const expanded = new Set<string>()
  for (const [name, steps] of [...abridged].sort((a, b) => a[1] - b[1])) {
    if (kept + steps > stepBudget) {
      break
    }
    kept += steps
    expanded.add(name)
  }
  return expanded
}
