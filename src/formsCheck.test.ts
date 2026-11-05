import { readFileSync } from 'node:fs'

import { formContradictions } from './formsCheck.ts'
import {
  NO_CALL_STATE,
  RARER_STATE,
  structuralForms,
} from './structuralForms.ts'

import type { StructuralForms } from './structuralForms.ts'
import type { Departure, HaplotypeDepartures } from './walkDepartures.ts'

const WINDOW = { start: 0, end: 100_000 }

const dep = (
  refStart: number,
  refEnd: number,
  altLength: number,
  more: Partial<Departure> = {},
): Departure => ({
  refStart,
  refEnd,
  altLength,
  invertedBp: 0,
  repeatedBp: 0,
  ...more,
})

const hap = (
  haplotype: string,
  departures: Departure[] = [],
  covered: [number, number][] = [[WINDOW.start, WINDOW.end]],
): HaplotypeDepartures => ({ haplotype, covered, departures })

const haps = (prefix: string, count: number, departures: Departure[] = []) =>
  Array.from({ length: count }, (_, i) => hap(`${prefix}${i}#1`, departures))

const kinds = (input: HaplotypeDepartures[], result: StructuralForms) =>
  formContradictions(input, result, { window: WINDOW }).map(c => [
    c.kind,
    c.haplotype,
  ])

// One haplotype moved from the form it was placed in to another
const moved = (result: StructuralForms, name: string, to: number) => ({
  ...result,
  forms: result.forms.map((f, i) => ({
    ...f,
    members:
      i === to ? [...f.members, name] : f.members.filter(m => m !== name),
  })),
})

const input = [
  ...haps('R', 6),
  ...haps('D', 6, [dep(10_000, 30_000, 0)]),
  ...haps('I', 6, [dep(12_000, 12_000, 14_000)]),
  ...haps('B', 3, [dep(25_000, 25_000, 14_000)]),
  hap('X#1', [dep(10_000, 30_000, 0), dep(60_000, 60_000, 400)]),
  hap(
    'G#1',
    [],
    [
      [0, 11_000],
      [29_000, 100_000],
    ],
  ),
  ...haps('V', 3, [dep(70_000, 70_300, 0), dep(70_400, 70_400, 260)]),
  ...haps('W', 3, [dep(70_000, 70_500, 0)]),
]

test('structuralForms states nothing its members contradict, at any floor or fold', () => {
  for (const floorBp of [50, 1000, 10_000]) {
    for (const foldRareForms of [true, false]) {
      for (const minCarriers of [1, 3, 5]) {
        const result = structuralForms(input, {
          window: WINDOW,
          floorBp,
          foldRareForms,
          minCarriers,
        })
        expect(kinds(input, result)).toEqual([])
      }
    }
  }
})

test('a member in a form whose change it lacks contradicts it in count', () => {
  const result = structuralForms(input, {
    window: WINDOW,
    floorBp: 50,
    minCarriers: 3,
  })
  const deletion = result.forms.findIndex(f => f.members.includes('D0#1'))
  expect(kinds(input, moved(result, 'R0#1', deletion))).toEqual([
    ['count', 'R0#1'],
  ])
})

test('a member shown with a call where its walks leave a gap contradicts it', () => {
  const result = structuralForms(input, {
    window: WINDOW,
    floorBp: 50,
    minCarriers: 3,
  })
  const reference = result.forms.findIndex(f => f.members.includes('R0#1'))
  expect(kinds(input, moved(result, 'G#1', reference))).toEqual([
    ['call', 'G#1'],
  ])
})

test('a member whose change lies elsewhere in the site contradicts the place', () => {
  const result = structuralForms(input, {
    window: WINDOW,
    floorBp: 50,
    minCarriers: 3,
  })
  const inserted = result.forms.findIndex(f => f.members.includes('I0#1'))
  expect(kinds(input, moved(result, 'B0#1', inserted))).toEqual([
    ['place', 'B0#1'],
  ])
})

test('a summed site states a net change, so pieces that cancel agree with ref', () => {
  const result = structuralForms(input, {
    window: WINDOW,
    floorBp: 50,
    minCarriers: 3,
  })
  expect(result.forms.find(f => f.members.includes('V0#1'))!.key).toBe(
    result.forms.find(f => f.members.includes('R0#1'))!.key,
  )
  expect(kinds(input, result)).toEqual([])
})

test("a shift's net agrees to the last figure its key gives it", () => {
  const shift = (net: number) => [
    dep(1000, 1000, 32_858),
    dep(12_000, 44_737 + net, 0),
  ]
  const shifted = [
    ...haps('R', 4),
    ...haps('S', 3, shift(0)),
    ...haps('T', 2, shift(-120)),
  ]
  const result = structuralForms(shifted, {
    window: WINDOW,
    floorBp: 50,
    minCarriers: 3,
  })
  expect(result.forms.map(f => f.members.length)).toEqual([5, 4])
  expect(kinds(shifted, result)).toEqual([])
})

// A result to hand-edit into a wrong one: `form` finds a member's form and
// `site` a site by where it starts
const edited = () => {
  const result = structuredClone(
    structuralForms(input, { window: WINDOW, floorBp: 50, minCarriers: 3 }),
  )
  const form = (name: string) =>
    result.forms.find(f => f.members.includes(name))!
  const site = (start: number) =>
    result.sites.findIndex(s => s.refStart === start)
  const state = (start: number, name: string) => {
    const i = site(start)
    const k = result.informative.indexOf(i)
    return result.sites[i]!.states.find(s => s.key === form(name).states[k])!
  }
  return { result, form, site, state }
}

const reported = (result: StructuralForms) => [
  ...new Set(kinds(input, result).map(([kind]) => kind)),
]

test('the checker reports a no call stated where the walks span the site', () => {
  const { result, form, site } = edited()
  form('D0#1').states[result.informative.indexOf(site(10_000))] = NO_CALL_STATE
  expect(reported(result)).toEqual(['hidden'])
})

test('the checker reports a rarer change stated where the member makes none', () => {
  const { result, form, site } = edited()
  form('R0#1').states[result.informative.indexOf(site(10_000))] = RARER_STATE
  expect(reported(result)).toEqual(['rarer'])
})

test("the checker reports a summed site's net change the members do not make", () => {
  const { result, state } = edited()
  state(70_000, 'W0#1').marks[0]!.delta = -1500
  expect(reported(result)).toEqual(['net'])
})

test("the checker reports a summed site's largest change the members do not make", () => {
  const { result, state } = edited()
  const [m] = state(70_000, 'W0#1').marks
  state(70_000, 'W0#1').marks = [
    { ...m!, kind: 'deletion', delta: -1500 },
    { ...m!, kind: 'insertion', refEnd: m!.refStart, delta: 1000 },
  ]
  expect(reported(result)).toEqual(['size'])
})

test("the checker reports a wider site's change stated at its place with another net change", () => {
  const { result, state } = edited()
  state(10_000, 'D0#1').marks[0]!.delta = -30_000
  expect(reported(result)).toEqual(['net'])
})

test("the checker reports a wider site's change stated at its place and net with another size", () => {
  const { result, state } = edited()
  const m = state(10_000, 'D0#1').marks[0]!
  state(10_000, 'D0#1').marks[0] = { ...m, kind: 'shift', shiftedBp: 40_000 }
  expect(reported(result)).toEqual(['size'])
})

test('the checker reports a change that lies in no site', () => {
  const { result, site } = edited()
  const i = site(60_000)
  result.sites[i] = { ...result.sites[i]!, refStart: 61_000, refEnd: 61_000 }
  expect(kinds(input, result)).toContainEqual(['unsited', 'X#1'])
})

test("the checker reports a member's changes elsewhere in a wider site than the form's, though they touch", () => {
  const shifted = [
    ...haps('R', 80),
    ...haps('A', 10, [dep(10_000, 20_000, 0)]),
    ...haps('B', 10, [dep(15_500, 25_500, 0)]),
  ]
  const result = structuralForms(shifted, { window: WINDOW, floorBp: 50 })
  expect(result.forms.map(f => f.key).sort()).toEqual([
    '-10000@10000',
    '-10000@15500',
    '0',
  ])
  expect(formContradictions(shifted, result, { window: WINDOW })).toEqual([])
  const into = result.forms.findIndex(f => f.members.includes('A0#1'))
  expect(
    formContradictions(shifted, moved(result, 'B0#1', into), {
      window: WINDOW,
    }).map(c => [c.kind, c.haplotype]),
  ).toEqual([['place', 'B0#1']])
})

test("the checker reports a member's changes in a wider site that the form does not state, away from its marks", () => {
  const extra = [
    ...haps('R', 10),
    ...haps('D', 10, [dep(10_000, 20_000, 0)]),
    ...haps('E', 10, [dep(15_000, 35_000, 0)]),
    hap('X#1', [
      dep(10_000, 20_000, 0),
      dep(30_000, 30_000, 800),
      dep(33_000, 33_800, 0),
    ]),
  ]
  const result = structuralForms(extra, {
    window: WINDOW,
    floorBp: 50,
    minCarriers: 5,
  })
  expect(result.sites.map(s => [s.refStart, s.refEnd])).toEqual([
    [10_000, 35_000],
  ])
  const into = result.forms.findIndex(f => f.members.includes('D0#1'))
  expect(
    formContradictions(extra, moved(result, 'X#1', into), {
      window: WINDOW,
    }).map(c => [c.kind, c.haplotype]),
  ).toEqual([['count', 'X#1']])
})

// T#1 deletes 5 kb and inserts it again at the common 100 bp insertion's
// site: a shift, whose insertion is no change of its own there
test("a rare shift whose insertion lies at a summed site is that site's change, kept apart", () => {
  const shifted = [
    ...haps('R', 10),
    ...haps('S', 5, [dep(50_000, 50_000, 100)]),
    hap('T#1', [dep(40_000, 45_000, 0), dep(50_000, 50_000, 5000)]),
  ]
  const result = structuralForms(shifted, {
    window: WINDOW,
    floorBp: 50,
    minCarriers: 3,
    foldRareForms: false,
  })
  expect(
    result.forms.find(f => f.members.includes('T#1'))!.marks.map(m => m.kind),
  ).toEqual(['shift'])
  expect(formContradictions(shifted, result, { window: WINDOW })).toEqual([])
})

test('the checker takes only the tolerance and constants from structuralForms', () => {
  const source = readFileSync(new URL('formsCheck.ts', import.meta.url), 'utf8')
  const imported = [
    ...source.matchAll(/^import \{([^}]*)\} from '\.\/structuralForms\.ts'/gm),
  ].flatMap(([, names]) => names!.split(',').map(name => name.trim()))
  expect(
    imported.filter(name => name !== '' && name !== name.toUpperCase()),
  ).toEqual(['agrees'])
})
