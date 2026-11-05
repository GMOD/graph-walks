import { formContradictions } from './formsCheck.ts'
import { structuralForms } from './structuralForms.ts'

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

test('a member in a form whose change it lacks contradicts it in net', () => {
  const result = structuralForms(input, {
    window: WINDOW,
    floorBp: 50,
    minCarriers: 3,
  })
  const deletion = result.forms.findIndex(f => f.members.includes('D0#1'))
  expect(kinds(input, moved(result, 'R0#1', deletion))).toEqual([
    ['net', 'R0#1'],
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
