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
// R's bases lo and hi
type Event =
  | { kind: 'on'; lo: number; hi: number; reversed: boolean }
  | { kind: 'own'; bp: number }
  | { kind: 'abridged'; lo: number; hi: number; bp: number }

// A point on a path's contig: where it is on R until (a fragment's last step
// or a departure's leave, `order` 0), from (a fragment's first step or a
// departure's rejoin, 1), or where an abridged stretch starts (2). `own` is
// the bp of the path's own that follow a leave, or precede a rejoin whose
// leave lies in an abridged stretch
interface Anchor {
  hap: number
  order: 0 | 1 | 2
  pos: number
  reversed: boolean
  exit?: true
  own?: number
  abridged?: { lo: number; hi: number; bp: number }
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
 * no call. An abridged fragment's bp, with those of the fragments and
 * departures lying within them, are a net change between the R its handles
 * reach. The walks on R's own haplotype are left out.
 */
export function walkDepartures(
  { fragments, departures, nodes }: DepartureSource,
  refName: string,
) {
  const refStart = (id: number) => {
    const node = nodes.get(id)
    return node?.rank === 0 && node.refName === refName ? node.start : NaN
  }
  const length = (id: number) => {
    const node = nodes.get(id)
    return node ? node.end - node.start : 0
  }
  const reference = panSNHaplotype(refName) ?? refName
  const byPath = new Map<string, Item[]>()
  const add = (name: string, item: Item) => {
    byPath.set(name, [...(byPath.get(name) ?? []), item])
  }
  for (const f of fragments) {
    add(f.name, fragmentItem(f, refStart, length))
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
    for (const chain of pathChains(items, refStart, length)) {
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

function fragmentItem(
  f: WalkFragment,
  refStart: (id: number) => number,
  length: (id: number) => number,
): Item {
  const n = f.ids.length
  const a = f.ids[0]!
  const b = f.ids[n - 1]!
  const sa = refStart(a)
  const sb = refStart(b)
  const lo = least(sa, sb)
  const hi = most(sa + length(a), sb + length(b))
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
 * departures leaving its last step and rejoining its first, fragments that
 * abut, and a departure touching no fragment on its own. What the read did
 * not see lies between chains, so each reads on its own, and R that two
 * chains reach has no call. Abridged fragments, each a net over the pieces
 * a writer moved out of a chunk, take in whatever of the path lies within
 * those pieces: the fragments of pieces in other chunks between them,
 * abridged fragments abutting them, and the departures whose both ends lie
 * in them, whose own steps the pieces hold and whose missing bp they do not.
 * A departure with one end in them keeps the bp of its own outside: those
 * the graph lacks and the steps the cut holds.
 */
function pathChains(
  items: Item[],
  refStart: (id: number) => number,
  length: (id: number) => number,
) {
  const overlaps = (a: Item, b: Item) =>
    a.hap[0] < b.hap[1] && a.hap[1] > b.hap[0]
  const touches = (a: Item, b: Item) =>
    a.hap[0] <= b.hap[1] && a.hap[1] >= b.hap[0]
  const extend = (g: Item, item: Item, bp: number) => {
    g.hap = [Math.min(g.hap[0], item.hap[0]), Math.max(g.hap[1], item.hap[1])]
    g.lo = least(g.lo, item.lo)
    g.hi = most(g.hi, item.hi)
    g.bp += bp
  }
  const abridged = (item: Item) =>
    'fragment' in item && item.fragment.abridged !== undefined
  let groups = items.filter(abridged).map(i => ({ ...i }))
  let rest = items.filter(i => !abridged(i))
  // the group whose pieces hold the step at `hap` on the contig, or whose
  // first or last step it is: own bp beside an abridged stretch join its net
  const holder = (hap: number) =>
    groups.find(g => g.hap[0] <= hap && hap <= g.hap[1])
  for (let changed = groups.length > 0; changed;) {
    const before = rest.length + groups.length
    rest = rest.filter(item => {
      if ('fragment' in item) {
        const g = groups.find(g => overlaps(item, g))
        if (g) {
          extend(g, item, item.bp)
        }
        return !g
      }
      const d = item.departure
      const a = holder(d.hapStart)
      const b = holder(d.hapStart + d.ownBp)
      if (a && b) {
        if (b !== a) {
          extend(a, b, b.bp)
          groups = groups.filter(g => g !== b)
        }
        extend(a, item, d.missingBp)
      }
      return !(a && b)
    })
    const merged: Item[] = []
    for (const g of groups.toSorted((a, b) => a.hap[0] - b.hap[0])) {
      const last = merged.at(-1)
      if (last && touches(last, g)) {
        extend(last, g, g.bp)
      } else {
        merged.push(g)
      }
    }
    groups = merged
    changed = rest.length + groups.length < before
  }
  // where each item starts and ends on the contig
  const edges = (item: Item): [number, number] =>
    'fragment' in item && !item.fragment.abridged
      ? [item.fragment.hapStart, item.fragment.hapEnd]
      : item.hap
  const all = [...groups, ...rest]
  const parent = new Map(all.map(item => [item, item]))
  const find = (item: Item): Item => {
    const p = parent.get(item) ?? item
    return p === item ? item : find(p)
  }
  const union = (a: Item, b: Item) => parent.set(find(a), find(b))
  // items meet where one's contig stretch touches another's: a departure's
  // leave at a fragment's last step or inside it, its rejoin at the first,
  // or two fragments abutting
  const meets = (a: [number, number], b: [number, number]) =>
    a[0] <= b[1] && a[1] >= b[0]
  all.forEach((item, i) => {
    for (const other of all.slice(i + 1)) {
      if (meets(edges(item), edges(other))) {
        union(item, other)
      }
    }
  })
  const chains = new Map<Item, Item[]>()
  for (const item of all) {
    const root = find(item)
    chains.set(root, [...(chains.get(root) ?? []), item])
  }
  return [...chains.values()].map(chain =>
    chainEvents(chain, holder, refStart, length),
  )
}

function chainEvents(
  items: Item[],
  holder: (hap: number) => Item | undefined,
  refStart: (id: number) => number,
  length: (id: number) => number,
) {
  const anchors: Anchor[] = []
  for (const item of items) {
    if ('fragment' in item) {
      const f = item.fragment
      if (f.abridged) {
        if (onR(item.lo) && onR(item.hi)) {
          anchors.push({
            hap: item.hap[0],
            order: 2,
            pos: item.lo,
            reversed: false,
            abridged: { lo: item.lo, hi: item.hi, bp: item.bp },
          })
        }
        continue
      }
      const n = f.ids.length
      const a = f.ids[0]!
      const b = f.ids[n - 1]!
      const ra = f.rev[0] === 1
      const rb = f.rev[n - 1] === 1
      if (onR(refStart(a))) {
        const pos = ra ? refStart(a) + length(a) : refStart(a)
        anchors.push({ hap: f.hapStart, order: 1, pos, reversed: ra })
      }
      if (onR(refStart(b))) {
        const pos = rb ? refStart(b) : refStart(b) + length(b)
        anchors.push({ hap: f.hapEnd, order: 0, pos, reversed: rb, exit: true })
      }
    } else {
      const d = item.departure
      const leaveIn = holder(d.hapStart) !== undefined
      const rejoinIn = holder(d.hapStart + d.ownBp) !== undefined
      const own =
        leaveIn || rejoinIn
          ? d.missingBp +
            (d.steps?.ids.reduce((sum, id) => sum + length(id), 0) ?? 0)
          : d.ownBp
      if (!leaveIn) {
        anchors.push({
          hap: d.hapStart,
          order: 0,
          pos: d.leavePos,
          reversed: d.leaveReversed,
          own,
        })
      }
      if (!rejoinIn) {
        anchors.push({
          hap: d.hapStart + d.ownBp,
          order: 1,
          pos: d.rejoinPos,
          reversed: d.rejoinReversed,
          ...(leaveIn ? { own } : {}),
        })
      }
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
      if (a.own !== undefined) {
        events.push({ kind: 'own', bp: a.own })
      }
      cur ??= { pos: a.pos, reversed: a.reversed }
    } else if (a.exit) {
      end = a.pos
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
    if (e.kind === 'on') {
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
  // its handles reach
  let abridgedEnd = NaN
  for (const e of flip ? events.toReversed() : events) {
    if (e.kind === 'own') {
      alt += e.bp
    } else if (e.kind === 'abridged') {
      if (Number.isNaN(prevEnd)) {
        first = e.lo
        prevEnd = e.lo
      }
      alt += e.bp
      abridgedEnd = e.hi
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
          invertedBp: 0,
          repeatedBp: 0,
          abridged: true,
        })
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
        })
      }
      prevEnd = e.hi
      reach = Math.max(reach, e.hi)
      alt = 0
      inverted = 0
    }
  }
  if (!Number.isNaN(abridgedEnd)) {
    h.departures.push({
      refStart: prevEnd,
      refEnd: Math.max(abridgedEnd, prevEnd),
      altLength: alt,
      invertedBp: 0,
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
  stepBudget?: number
  // the most compressed bytes either cut fetches from the region's chunks,
  // past which the cut throws RegionTooLargeError; no limit by default. Only
  // stepBudget bounds the detail records expansion reads
  byteBudget?: number
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
    signal,
  } = opts
  const header: WalkReaderHeader = await reader.header({ signal })
  const sample = header.references[0] ?? ''
  const refName =
    (await reader.sequenceFor(sample, region.refName, { signal })) ??
    region.refName
  const cutOpts = { links: false, byteBudget, signal }
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
