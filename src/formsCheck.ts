import {
  NO_CALL_STATE,
  OTHER_STATE,
  RARER_STATE,
  REFERENCE_STATE,
  SHIFT_MIN_BP,
  SHIFT_TOLERANCE,
  STRUCTURAL_BP,
  SUMMED_SITE_BP,
  agrees,
} from './structuralForms.ts'

import type { Form, Site, StructuralForms } from './structuralForms.ts'
import type { Departure, HaplotypeDepartures } from './walkDepartures.ts'

/**
 * How a form's statement at a site can contradict a member's own walks:
 * - `call`: the form states a call where the member's walks do not span the
 *   site's part inside the window
 * - `hidden`: the form states no call where they do
 * - `rarer`: the form states a rarer change where the member makes none
 * - `count`: at a site wider than SUMMED_SITE_BP, the member makes more or
 *   fewer changes there than the form states
 * - `place`: at such a site, as many, but they and the form's marks pair
 *   one-to-one on no pairing of changes that cover half of each other: one
 *   lies elsewhere in the site
 * - `net`: at a summed site, the member's net change there and the form's
 *   disagree; at a wider one, such pairings exist, but none pairs only
 *   changes agreeing in net change
 * - `size`: the same for the largest change at a summed site, and for each
 *   pair's size at a wider one
 * - `unsited`: a change of the member's lies in no site
 */
export type ContradictionKind =
  'call' | 'net' | 'size' | 'rarer' | 'count' | 'place' | 'unsited' | 'hidden'

export interface Contradiction {
  kind: ContradictionKind
  haplotype: string
  site: number
  stated: number
  own: number
}

// A run of a haplotype's departures, joined where nearer than joinBp on R
interface Piece {
  refStart: number
  refEnd: number
  delta: number
  invertedBp: number
  repeatedBp: number
  abridged: boolean
}

// A change a haplotype makes, as the tolerance compares it: its net change
// and its size, the span of an inversion and the deletion of a shift.
// `plain` changes insert, delete or replace
interface Change {
  refStart: number
  refEnd: number
  delta: number
  size: number
  plain: boolean
}

const touches = (
  m: { refStart: number; refEnd: number },
  start: number,
  end: number,
) => m.refStart <= end && m.refEnd >= start

// Whether each of two changes covers at least half of the other on R, a
// point taking 1 bp
const halfCover = (a: Change, b: Change) => {
  const aEnd = Math.max(a.refEnd, a.refStart + 1)
  const bEnd = Math.max(b.refEnd, b.refStart + 1)
  const overlap = Math.min(aEnd, bEnd) - Math.max(a.refStart, b.refStart)
  return 2 * overlap >= aEnd - a.refStart && 2 * overlap >= bEnd - b.refStart
}

const sum = (ns: number[]) => ns.reduce((s, n) => s + n, 0)

const largest = (changes: Change[]) =>
  changes.reduce((most, c) => Math.max(most, c.size), 0)

function joined(departures: Departure[], joinBp: number) {
  const out: Piece[] = []
  for (const d of departures.toSorted(
    (a, b) => a.refStart - b.refStart || a.refEnd - b.refEnd,
  )) {
    const delta = d.altLength - (d.refEnd - d.refStart)
    const last = out.at(-1)
    if (
      last &&
      !last.abridged &&
      !d.abridged &&
      d.refStart - last.refEnd < joinBp
    ) {
      last.refEnd = Math.max(last.refEnd, d.refEnd)
      last.delta += delta
      last.invertedBp += d.invertedBp
      last.repeatedBp += d.repeatedBp
    } else {
      out.push({
        refStart: d.refStart,
        refEnd: d.refEnd,
        delta,
        invertedBp: d.invertedBp,
        repeatedBp: d.repeatedBp,
        abridged: d.abridged === true,
      })
    }
  }
  return out
}

/**
 * A haplotype's changes of at least `floorBp` over the window, from its
 * departures: pieces joined where nearer than joinBp; one changing R's
 * length by floorBp or more changes it, an abridged one by its net, and one
 * that does not but walks floorBp of R against it inverts that span. A
 * deletion of SHIFT_MIN_BP or more pairs with the first insertion or step
 * back in R order within SHIFT_TOLERANCE of its size and twice its size of
 * it, largest deletion first, into one shift. `plain` lists the pieces a
 * summed site adds up: those neither abridged nor making a change other than
 * an unpaired insertion, deletion or replacement.
 */
function haplotypeChanges(
  departures: Departure[],
  floorBp: number,
  joinBp: number,
  window: { start: number; end: number },
) {
  const pieces = joined(departures, joinBp).filter(p =>
    touches(p, window.start, window.end),
  )
  const kinds = pieces.map(p => {
    const span = p.refEnd - p.refStart
    const alt = Math.max(0, span + p.delta)
    const size = Math.abs(p.delta)
    return p.abridged
      ? size >= floorBp
        ? 'net'
        : undefined
      : size < floorBp
        ? p.invertedBp >= floorBp
          ? 'inversion'
          : undefined
        : p.repeatedBp >= STRUCTURAL_BP && p.delta > 0
          ? 'duplication'
          : alt > span && span < STRUCTURAL_BP
            ? 'insertion'
            : alt < span && alt < STRUCTURAL_BP
              ? 'deletion'
              : 'replacement'
  })
  const paired = new Map<number, number>()
  const deletions = pieces
    .map((_, i) => i)
    .filter(i => kinds[i] === 'deletion' && -pieces[i]!.delta >= SHIFT_MIN_BP)
    .sort((a, b) => pieces[a]!.delta - pieces[b]!.delta)
  for (const i of deletions) {
    const d = pieces[i]!
    const size = -d.delta
    const gain = pieces.findIndex(
      (g, j) =>
        !paired.has(j) &&
        (kinds[j] === 'insertion' || kinds[j] === 'duplication') &&
        Math.abs(g.delta - size) <= SHIFT_TOLERANCE * Math.max(g.delta, size) &&
        Math.max(g.refStart - d.refEnd, d.refStart - g.refEnd) <= 2 * size,
    )
    if (gain >= 0) {
      paired.set(gain, i)
      paired.set(i, gain)
    }
  }
  const changes: Change[] = []
  const plain: Piece[] = []
  pieces.forEach((p, i) => {
    const kind = kinds[i]
    const other = paired.get(i)
    const isPlain =
      kind === 'insertion' || kind === 'deletion' || kind === 'replacement'
    if (other !== undefined) {
      if (kind === 'deletion') {
        const g = pieces[other]!
        changes.push({
          refStart: Math.min(p.refStart, g.refStart),
          refEnd: Math.max(p.refEnd, g.refEnd),
          delta: p.delta + g.delta,
          size: -p.delta,
          plain: false,
        })
      }
    } else if (kind === undefined) {
      if (!p.abridged) {
        plain.push(p)
      }
    } else {
      changes.push({
        refStart: p.refStart,
        refEnd: p.refEnd,
        delta: p.delta,
        size: kind === 'inversion' ? p.refEnd - p.refStart : Math.abs(p.delta),
        plain: isPlain,
      })
      if (isPlain) {
        plain.push(p)
      }
    }
  })
  return { changes, plain }
}

// A state's marks as changes, sized as a member's are
const statedChanges = (site: Site, key: string): Change[] =>
  (site.states.find(s => s.key === key)?.marks ?? []).map(m => ({
    refStart: m.refStart,
    refEnd: m.refEnd,
    delta: m.delta,
    size:
      m.kind === 'inversion'
        ? m.refEnd - m.refStart
        : (m.shiftedBp ?? Math.abs(m.delta)),
    plain:
      m.kind === 'insertion' ||
      m.kind === 'deletion' ||
      m.kind === 'replacement',
  }))

// Whether `as` and `bs` pair one-to-one with `pairs` true of every pair
function pairOneToOne<T>(as: T[], bs: T[], pairs: (a: T, b: T) => boolean) {
  const holder = bs.map(() => -1)
  const place = (i: number, seen: Set<number>): boolean =>
    bs.some((b, j) => {
      let placed = false
      if (!seen.has(j) && pairs(as[i]!, b)) {
        seen.add(j)
        placed = holder[j]! < 0 || place(holder[j]!, seen)
        if (placed) {
          holder[j] = i
        }
      }
      return placed
    })
  return as.length === bs.length && as.every((_, i) => place(i, new Set()))
}

// The state key a form states at a site: its key at an informative site, and
// at a shared one the state of its marks there, R's where it has none; none
// at an edge site or for `other`
function statedKey(result: StructuralForms, form: Form, site: number) {
  const k = result.informative.indexOf(site)
  return form.key === OTHER_STATE || result.sites[site]!.edge
    ? undefined
    : k >= 0
      ? form.states[k]
      : result.sites[site]!.shared
        ? (form.marks.find(m => m.site === site)?.state ?? REFERENCE_STATE)
        : undefined
}

/**
 * Every place a structuralForms result states something of a member that the
 * member's own departures contradict, from the departures and the result's
 * statements alone: at each site that is not `edge`, of each form but
 * `other`, the member's coverage against the form's call, and its changes
 * touching the site against the state's marks there. At a site at most
 * SUMMED_SITE_BP wide its plain pieces there add up to one change where they
 * reach the floor, and its net and largest change agree with the state's;
 * at a wider site its changes pair one-to-one with the state's marks, each
 * pair covering half of each other and agreeing in net change and size.
 * Agreement is structuralForms' `agrees`, the one definition the two share.
 */
export function formContradictions(
  input: HaplotypeDepartures[],
  result: StructuralForms,
  {
    window,
    joinBp = STRUCTURAL_BP,
  }: { window: { start: number; end: number }; joinBp?: number },
) {
  const { floorBp, sites } = result
  const byName = new Map(input.map(h => [h.haplotype, h]))
  const out: Contradiction[] = []
  const covers = (h: HaplotypeDepartures, site: Site) => {
    const start = Math.max(site.refStart, window.start)
    const end = Math.min(site.refEnd, window.end)
    return h.covered.some(([a, b]) => a <= start && b >= end)
  }
  const netAgrees = (a: Change, b: Change) =>
    halfCover(a, b) &&
    agrees(a.delta, b.delta, floorBp, Math.max(a.size, b.size))
  const sizeAgrees = (a: Change, b: Change) =>
    netAgrees(a, b) && agrees(a.size, b.size, floorBp)
  for (const form of result.forms) {
    for (const name of form.members) {
      const h = byName.get(name)
      const { changes, plain } = haplotypeChanges(
        h?.departures ?? [],
        floorBp,
        joinBp,
        window,
      )
      const add = (
        kind: ContradictionKind,
        site: number,
        stated: number,
        own: number,
      ) => {
        out.push({ kind, haplotype: name, site, stated, own })
      }
      for (const c of changes) {
        if (!sites.some(s => touches(c, s.refStart, s.refEnd))) {
          add('unsited', -1, 0, c.delta)
        }
      }
      sites.forEach((site, i) => {
        const key = statedKey(result, form, i)
        if (key !== undefined) {
          const called = h !== undefined && covers(h, site)
          const summed = site.refEnd - site.refStart <= SUMMED_SITE_BP
          const at = changes.filter(c => touches(c, site.refStart, site.refEnd))
          const plainNet = sum(
            plain
              .filter(p => touches(p, site.refStart, site.refEnd))
              .map(p => p.delta),
          )
          const own = summed
            ? [
                ...at.filter(c => !c.plain),
                ...(Math.abs(plainNet) >= floorBp
                  ? [
                      {
                        refStart: site.refStart,
                        refEnd: site.refEnd,
                        delta: plainNet,
                        size: Math.abs(plainNet),
                        plain: true,
                      },
                    ]
                  : []),
              ]
            : at
          const ownNet = sum(own.map(c => c.delta))
          const stated = statedChanges(site, key)
          const statedNet = sum(stated.map(c => c.delta))
          if (key === NO_CALL_STATE) {
            if (called) {
              add('hidden', i, 0, ownNet)
            }
          } else if (!called) {
            add('call', i, 0, 0)
          } else if (key === RARER_STATE) {
            if (own.length === 0) {
              add('rarer', i, floorBp, ownNet)
            }
          } else if (summed) {
            if (
              !agrees(
                statedNet,
                ownNet,
                floorBp,
                Math.max(largest(stated), largest(own)),
              )
            ) {
              add('net', i, statedNet, ownNet)
            } else if (!agrees(largest(stated), largest(own), floorBp)) {
              add('size', i, largest(stated), largest(own))
            }
          } else if (stated.length !== own.length) {
            add('count', i, stated.length, own.length)
          } else if (!pairOneToOne(stated, own, halfCover)) {
            add('place', i, statedNet, ownNet)
          } else if (!pairOneToOne(stated, own, netAgrees)) {
            add('net', i, statedNet, ownNet)
          } else if (!pairOneToOne(stated, own, sizeAgrees)) {
            add('size', i, largest(stated), largest(own))
          }
        }
      })
    }
  }
  return out
}
