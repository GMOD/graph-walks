import {
  NO_CALL_STATE,
  OTHER_STATE,
  RARER_STATE,
  STRUCTURAL_BP,
  SUMMED_SITE_BP,
  agrees,
  formStateAt,
  haplotypeMarks,
  joinDepartures,
  largestChange,
  statesAgree,
} from './structuralForms.ts'

import type { Mark, Site, StructuralForms } from './structuralForms.ts'
import type { HaplotypeDepartures } from './walkDepartures.ts'

/**
 * How a form's statement at a site can contradict a member's own walks:
 * - `call`: the form states a call where the member's walks do not span the
 *   site's part inside the window
 * - `net`: the member's net change there differs from the form's
 * - `size`: its largest change there differs, its net standing for its plain
 *   changes at a site at most SUMMED_SITE_BP wide
 * - `rarer`: the form states a rarer change where the member makes none
 * - `place`: at a site wider than SUMMED_SITE_BP, the member's marks touching
 *   one of the form's there differ from the form's marks touching it
 * - `unsited`: a change of the member's lies in no site
 * - `hidden`: the form states no call where the member has one
 */
export type ContradictionKind =
  'call' | 'net' | 'size' | 'rarer' | 'place' | 'unsited' | 'hidden'

export interface Contradiction {
  kind: ContradictionKind
  haplotype: string
  site: number
  stated: number
  own: number
}

const touches = (
  m: { refStart: number; refEnd: number },
  start: number,
  end: number,
) => m.refStart <= end && m.refEnd >= start

const PLAIN_KINDS = new Set(['insertion', 'deletion', 'replacement'])

const sum = (ns: number[]) => ns.reduce((s, n) => s + n, 0)

/**
 * Every place a structuralForms result states something of a member that the
 * member's own departures contradict, read off the departures, not off the
 * sites' states: at each site that is not `edge`, the member's net change
 * (its departures, joined, touching the site where the site is at most
 * SUMMED_SITE_BP wide, its marks there otherwise) and its largest mark
 * against the form's marks touching the site, and its coverage against the
 * form's call. An `other` member's form states nothing.
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
  for (const form of result.forms) {
    for (const name of form.members) {
      const h = byName.get(name)!
      const departures = joinDepartures(h.departures, joinBp).filter(d =>
        touches(d, window.start, window.end),
      )
      const marks = haplotypeMarks(h.departures, { floorBp, joinBp, window })
      const add = (
        kind: ContradictionKind,
        site: number,
        stated: number,
        own: number,
      ) => {
        out.push({ kind, haplotype: name, site, stated, own })
      }
      for (const m of marks) {
        if (!sites.some(s => touches(m, s.refStart, s.refEnd))) {
          add('unsited', -1, 0, m.delta)
        }
      }
      sites.forEach((site, i) => {
        const state = site.edge ? OTHER_STATE : formStateAt(result, form, i)
        if (state !== OTHER_STATE && state !== undefined) {
          const called = covers(h, site)
          const at = <T extends { refStart: number; refEnd: number }>(
            xs: T[],
          ) => xs.filter(x => touches(x, site.refStart, site.refEnd))
          const ownMarks = at(marks)
          const summed = site.refEnd - site.refStart <= SUMMED_SITE_BP
          const ownNet = summed
            ? sum(at(departures).map(d => d.delta))
            : sum(ownMarks.map(m => m.delta))
          const ownSize = summed
            ? Math.max(
                Math.abs(ownNet),
                largestChange(ownMarks.filter(m => !PLAIN_KINDS.has(m.kind))),
              )
            : largestChange(ownMarks)
          if (state === NO_CALL_STATE) {
            if (called) {
              add('hidden', i, 0, ownNet)
            }
          } else if (!called) {
            add('call', i, 0, 0)
          } else if (state === RARER_STATE) {
            if (ownSize < floorBp) {
              add('rarer', i, floorBp, ownNet)
            }
          } else {
            const stated = site.states.find(st => st.key === state)?.marks ?? []
            const statedNet = sum(stated.map(m => m.delta))
            if (
              !agrees(
                statedNet,
                ownNet,
                floorBp,
                Math.max(largestChange(stated), ownSize),
              )
            ) {
              add('net', i, statedNet, ownNet)
            } else if (!agrees(largestChange(stated), ownSize, floorBp)) {
              add('size', i, largestChange(stated), ownSize)
            } else if (!summed) {
              for (const c of stated) {
                const near = <T extends Mark>(xs: T[]) =>
                  xs.filter(x => touches(x, c.refStart, c.refEnd))
                if (!statesAgree(near(stated), near(ownMarks), floorBp)) {
                  add(
                    'place',
                    i,
                    c.delta,
                    sum(near(ownMarks).map(m => m.delta)),
                  )
                }
              }
            }
          }
        }
      })
    }
  }
  return out
}
