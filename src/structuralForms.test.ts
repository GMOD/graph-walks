import {
  NO_CALL_STATE,
  RARER_STATE,
  REFERENCE_STATE,
  autoFloor,
  changeKind,
  describeForm,
  formStateAt,
  formTag,
  haplotypeMarks,
  markTag,
  minCarriersFor,
  pairShifts,
  stateKey,
  structuralForms,
} from './structuralForms.ts'

import type { Mark } from './structuralForms.ts'
import type { Departure, HaplotypeDepartures } from './walkDepartures.ts'

const WINDOW = { start: 0, end: 100_000 }

// A departure over R's start–end with `alt` bp of the haplotype's own
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
): HaplotypeDepartures => ({ haplotype, covered, departures, open: false })

// `count` haplotypes named from `prefix`, each with the same departures
const haps = (prefix: string, count: number, departures: Departure[] = []) =>
  Array.from({ length: count }, (_, i) => hap(`${prefix}${i}#1`, departures))

const mark = (
  kind: Mark['kind'],
  refStart: number,
  refEnd: number,
  delta: number,
): Mark => ({
  kind,
  refStart,
  refEnd,
  altLength: refEnd - refStart + delta,
  delta,
})

const forms = (...input: Parameters<typeof structuralForms>) =>
  structuralForms(input[0], { window: WINDOW, floorBp: 50, ...input[1] })

const memberCounts = (r: ReturnType<typeof structuralForms>) =>
  r.forms.map(f => f.members.length)

test('a change under 50 bp is no mark, a larger one keys on two figures', () => {
  expect(haplotypeMarks([dep(100, 146, 0)], {})).toEqual([])
  const key = (delta: number) => stateKey(mark('deletion', 0, 0, delta))
  expect(key(-1716)).toBe('-1700')
  expect(key(-1740)).toBe('-1700')
  expect(key(-84_684)).toBe('-85000')
  expect(key(5600)).toBe('+5600')
  expect(key(22_201)).toBe(key(22_204))
  expect(key(11_200)).not.toBe(key(5600))
})

test('an inversion moves no length and is a state of its own', () => {
  const marks = haplotypeMarks([dep(1000, 1300, 300, { invertedBp: 300 })], {})
  expect(marks).toEqual([mark('inversion', 1000, 1300, 0)])
  expect(stateKey(marks[0]!)).toBe('v300')
})

test('a form needs 1% of the haplotypes, and at least one', () => {
  expect(minCarriersFor(462)).toBe(5)
  expect(minCarriersFor(11)).toBe(1)
})

// Checked against the release's wave callset on 2026-10-21: GSTM1's deletion
// is REF 18,446 to 1, C4's 26 kb gain REF 32,738 to 59,108, and LPA's KIV-2
// gains replace a 30,736 bp array
test('a change inserts, deletes or replaces on its exact REF and ALT', () => {
  expect(changeKind(1, 18_001)).toBe('insertion')
  expect(changeKind(18_446, 1)).toBe('deletion')
  expect(changeKind(32_738, 59_108)).toBe('replacement')
  expect(changeKind(30_736, 108_736)).toBe('replacement')
  expect(changeKind(176_380, 82_380)).toBe('replacement')
  expect(changeKind(100, 200)).toBe('replacement')
  expect(changeKind(5, 305)).toBe('insertion')
  // what the sidecar's rounding slack called a deletion
  expect(changeKind(18_446, 60)).toBe('replacement')
})

test('departures under 50 bp of R apart join, keeping their net change', () => {
  expect(
    haplotypeMarks(
      [dep(1000, 1400, 0), dep(1410, 1410, 2), dep(1420, 1500, 0)],
      {},
    ),
  ).toEqual([
    {
      kind: 'deletion',
      refStart: 1000,
      refEnd: 1500,
      altLength: 22,
      delta: -478,
    },
  ])
  expect(
    haplotypeMarks([dep(1000, 1400, 0), dep(1460, 1800, 0)], {}),
  ).toHaveLength(2)
})

test('forms group the haplotypes that match at every informative site', () => {
  const del = dep(100, 1800, 0)
  const result = forms(
    [
      hap('A#1', [del, dep(5000, 6700, 0)]),
      hap('A#2', [del]),
      hap('B#1', [del]),
      hap('B#2'),
      hap('C#1'),
      hap('C#2'),
    ],
    { minCarriers: 3 },
  )
  expect(result.informative).toEqual([0])
  expect(result.forms.map(f => [f.key, f.members])).toEqual([
    ['0', ['B#2', 'C#1', 'C#2']],
    ['-1700', ['A#1', 'A#2', 'B#1']],
  ])
  expect(result.rareCarriers).toEqual(['A#1'])
})

test('a form marks its changes in R order and says them largest first', () => {
  const result = forms(
    [
      hap('A#1', [dep(100, 1800, 0)]),
      hap('A#2', [dep(100, 1800, 0), dep(3000, 3000, 6400)]),
      hap('B#1', [dep(3000, 3000, 6400)]),
      hap('B#2', [
        dep(3000, 3000, 6400),
        dep(4000, 4500, 500, { invertedBp: 500 }),
      ]),
      hap('C#1', [dep(4000, 4500, 500, { invertedBp: 500 })]),
      hap('C#2'),
    ],
    { minCarriers: 1 },
  )
  const a2 = result.forms.find(f => f.members.includes('A#2'))!
  expect(a2.marks).toEqual([
    { ...mark('deletion', 100, 1800, -1700), site: 0, state: '-1700' },
    { ...mark('insertion', 3000, 3000, 6400), site: 1, state: '+6400' },
  ])
  expect(describeForm(a2)).toBe('6.4 kb insertion, 1.7 kb deletion')
  expect(describeForm(result.forms.find(f => f.members.includes('B#2'))!)).toBe(
    '6.4 kb insertion, 500 bp inversion',
  )
  expect(describeForm(result.forms.find(f => f.members.includes('C#2'))!)).toBe(
    'as the reference',
  )
})

test('more than two size changes are counted, the largest named', () => {
  expect(
    describeForm({
      states: [],
      marks: [
        mark('deletion', 0, 1000, -1000),
        mark('insertion', 2000, 2000, 19_000),
        mark('deletion', 3000, 5000, -2000),
      ],
    }),
  ).toBe('3 size changes, largest a 19 kb insertion')
})

// GSTM1: 77 haplotypes in the sidecar carry an 18 kb deletion over the gene
// and an 18 kb insertion at its end, which is the reference's length
test('a deletion and an insertion of one size pair into a shift', () => {
  expect(
    pairShifts([
      mark('deletion', 100, 18_200, -18_100),
      mark('insertion', 20_000, 20_000, 17_000),
      mark('insertion', 30_000, 30_000, 300),
    ]),
  ).toEqual([
    {
      kind: 'shift',
      refStart: 100,
      refEnd: 20_000,
      altLength: 19_900 - 1100,
      delta: -1100,
      shiftedBp: 18_100,
    },
    mark('insertion', 30_000, 30_000, 300),
  ])
  // over 10% apart in size, far apart on R, or a replacement: no shift
  expect(
    pairShifts([
      mark('deletion', 100, 18_200, -18_100),
      mark('insertion', 20_000, 20_000, 6400),
    ]),
  ).toHaveLength(2)
  expect(
    pairShifts([
      mark('deletion', 100, 18_200, -18_100),
      mark('insertion', 220_000, 220_000, 18_100),
    ]),
  ).toHaveLength(2)
  expect(
    pairShifts([
      mark('replacement', 0, 800_000, -8000),
      mark('insertion', 9000, 9000, 8000),
    ]).map(m => m.kind),
  ).toEqual(['replacement', 'insertion'])
})

test('a step back along R is a duplication, which shifts like an insertion', () => {
  const copy = dep(50_000, 50_000, 49_000, {
    repeatedBp: 49_000,
    copied: [1000, 50_000],
  })
  expect(haplotypeMarks([copy], {})).toEqual([
    {
      kind: 'duplication',
      refStart: 50_000,
      refEnd: 50_000,
      altLength: 49_000,
      delta: 49_000,
      repeatedBp: 49_000,
      copied: [1000, 50_000],
    },
  ])
  expect(describeForm({ states: [], marks: haplotypeMarks([copy], {}) })).toBe(
    '49 kb tandem duplication',
  )
  expect(
    haplotypeMarks([copy, dep(60_000, 108_000, 0)], {}).map(m => m.kind),
  ).toEqual(['shift'])
})

test('an abridged departure is a net mark, and its haplotype is listed', () => {
  const result = forms([
    hap('A#1', [dep(10_000, 75_000, 60_000, { abridged: true })]),
    hap('B#1'),
  ])
  expect(result.forms[1]!.marks).toEqual([
    { ...mark('net', 10_000, 75_000, -5000), site: 0, state: 'n-5000' },
  ])
  expect(describeForm(result.forms[1]!)).toBe('5 kb shorter in 65 kb')
  expect(result.abridged).toEqual(['A#1'])
})

test('a gain at a site wider than 50 bp replaces it', () => {
  const result = forms([
    ...haps('L', 6, [dep(1000, 31_000, 108_000)]),
    ...haps('S', 6, [dep(1000, 31_000, 10_000)]),
    ...haps('R', 8),
  ])
  expect(result.forms.map(f => describeForm(f))).toEqual([
    'as the reference',
    '30 kb replaced by 10 kb',
    '30 kb replaced by 108 kb',
  ])
})

test('a window with no change is one form', () => {
  const result = forms([hap('A#1'), hap('A#2')])
  expect(result.sites).toEqual([])
  expect(result.forms).toEqual([
    { key: '', states: [], members: ['A#1', 'A#2'], marks: [], net: 0 },
  ])
})

test('a change every haplotype shares tells none apart, and every form carries it', () => {
  const result = forms(haps('A', 4, [dep(100, 1800, 0)]), { minCarriers: 3 })
  expect(result.informative).toEqual([])
  expect(result.forms).toHaveLength(1)
  expect(describeForm(result.forms[0]!)).toBe('1.7 kb deletion')
})

test('nested changes are one site, each haplotype stating only its own', () => {
  const result = forms([
    ...haps('P', 3, [dep(1000, 21_000, 0)]),
    ...haps('C', 3, [dep(5000, 6000, 0)]),
    ...haps('R', 3),
  ])
  expect(result.sites.map(s => [s.refStart, s.refEnd])).toEqual([
    [1000, 21_000],
  ])
  expect(result.forms.map(f => [f.key, f.net])).toEqual([
    ['0', 0],
    ['-1000', -1000],
    ['-20000', -20_000],
  ])
})

test('a rare state folds to the common state nearest in net change', () => {
  const result = forms(
    [
      ...haps('R', 10),
      ...haps('D', 5, [dep(1000, 21_000, 0)]),
      hap('X#1', [dep(1000, 19_800, 0)]),
    ],
    { minCarriers: 3 },
  )
  expect(result.forms.map(f => [f.key, f.members.length])).toEqual([
    ['0', 10],
    ['-20000', 6],
  ])
  expect(result.rareCarriers).toEqual(['X#1'])
})

test('a form under minCarriers folds into the form nearest in net', () => {
  const big = dep(1000, 21_000, 0)
  const small = dep(40_000, 41_000, 0)
  const result = forms(
    [
      ...haps('R', 6),
      ...haps('B', 3, [big]),
      ...haps('S', 3, [small]),
      hap('X#1', [big, small]),
    ],
    { minCarriers: 2 },
  )
  expect(memberCounts(result)).toEqual([6, 4, 3])
  expect(result.placed).toBe(13)
  expect(result.forms[1]!.members).toContain('X#1')
  expect(result.rareCarriers).toEqual(['X#1'])
})

test('foldRareForms: false keeps a form under minCarriers apart', () => {
  const big = dep(1000, 21_000, 0)
  const small = dep(40_000, 41_000, 0)
  const result = forms(
    [
      ...haps('R', 6),
      ...haps('B', 3, [big]),
      ...haps('S', 3, [small]),
      hap('X#1', [big, small]),
    ],
    { minCarriers: 2, foldRareForms: false },
  )
  expect(memberCounts(result)).toEqual([6, 3, 3, 1])
  expect(result.rareCarriers).toEqual([])
})

test('a rare call where most haplotypes have no call stays apart', () => {
  const flank = dep(100, 400, 0)
  const del = (bp: number) => [dep(1000, 1000 + bp, 0)]
  const result = forms(
    [
      ...Array.from({ length: 5 }, (_, i) => hap(`N${i}#1`, [], [[0, 500]])),
      ...haps('R', 3, [flank]),
      hap('X#1', del(310)),
      hap('Y#1', del(620)),
    ],
    { minCarriers: 2 },
  )
  expect(result.forms.map(f => [f.key, f.members.length])).toEqual([
    ['0|.', 5],
    ['-300|0', 3],
    ['0|~', 2],
  ])
  expect(result.forms[2]!.states).toEqual([REFERENCE_STATE, RARER_STATE])
  expect(describeForm(result.forms[2]!)).toBe('a rarer change at 1 site')
  expect(describeForm(result.forms[0]!)).toBe('not aligned at 1 site')
})

test('a mark at a site the walks do not span whole is no call there', () => {
  const insertion = dep(50_000, 50_000, 5000)
  const result = forms(
    [
      ...haps('R', 4),
      ...haps('D', 3, [dep(10_000, 12_000, 0)]),
      ...haps('I', 4, [insertion]),
      ...haps('W', 3, [dep(49_000, 52_000, 9000)]),
      hap(
        'X#1',
        [insertion],
        [
          [0, 49_500],
          [50_000, 100_000],
        ],
      ),
    ],
    { minCarriers: 1 },
  )
  expect(result.sites.map(s => [s.refStart, s.refEnd])).toEqual([
    [10_000, 12_000],
    [49_000, 52_000],
  ])
  const x = result.forms.find(f => f.members.includes('X#1'))!
  expect(x.states).toEqual([REFERENCE_STATE, NO_CALL_STATE])
  expect(x.marks).toEqual([])
})

// AMY1's 12 in the sidecar: no call at the array, called in the flank
test('a haplotype placed somewhere is in a form, one placed nowhere is not', () => {
  const result = forms(
    [
      ...haps('R', 4, [dep(1000, 3000, 0)]),
      ...haps('D', 4, [dep(1000, 3000, 0), dep(60_000, 70_000, 0)]),
      ...Array.from({ length: 2 }, (_, i) =>
        hap(`L${i}#1`, [dep(1000, 3000, 0)], [[0, 30_000]]),
      ),
      hap('Z#1', [], [[200_000, 300_000]]),
    ],
    { minCarriers: 2 },
  )
  expect(result.forms.map(f => [f.key, f.members.length])).toEqual([
    ['0', 4],
    ['-10000', 4],
    ['.', 2],
  ])
  expect(describeForm(result.forms[2]!)).toBe(
    '2 kb deletion, not aligned at 1 site',
  )
  expect(result.placed).toBe(10)
  expect(result.unplaced).toEqual(['Z#1'])
})

test('each mark names its site and state, and formStateAt reads a form at any site', () => {
  const result = forms(
    [
      ...haps('R', 4, [dep(1000, 3000, 0)]),
      ...haps('D', 4, [dep(1000, 3000, 0), dep(60_000, 70_000, 0)]),
      ...Array.from({ length: 2 }, (_, i) =>
        hap(`L${i}#1`, [dep(1000, 3000, 0)], [[0, 30_000]]),
      ),
    ],
    { minCarriers: 2 },
  )
  expect(result.sites.map(s => [s.edge, s.shared])).toEqual([
    [false, true],
    [false, false],
  ])
  expect(result.informative).toEqual([1])
  expect(result.forms[1]!.marks.map(m => [m.site, m.state])).toEqual([
    [0, '-2000'],
    [1, '-10000'],
  ])
  expect(
    result.forms.map(f =>
      result.sites.map((_, s) => formStateAt(result, f, s)),
    ),
  ).toEqual([
    ['-2000', REFERENCE_STATE],
    ['-2000', '-10000'],
    ['-2000', NO_CALL_STATE],
  ])
  const edged = forms(
    [
      ...Array.from({ length: 6 }, (_, i) => hap(`S${i}#1`, [], [[0, 90_000]])),
      ...haps('E', 3, [dep(95_000, 99_000, 0)]),
    ],
    { minCarriers: 2 },
  )
  expect(edged.sites[0]!.shared).toBe(false)
  expect(formStateAt(edged, edged.forms[0]!, 0)).toBeUndefined()
})

test('markTag says a change in a few characters', () => {
  expect(markTag(mark('insertion', 0, 0, 78_000))).toBe('+78 kb')
  expect(markTag(mark('deletion', 0, 18_445, -18_445))).toBe('18 kb del')
  expect(markTag(mark('net', 0, 160_000, -94_000))).toBe('\u221294 kb')
  expect(markTag(mark('replacement', 0, 25_900, 78_100))).toBe('+78 kb')
  expect(markTag(mark('replacement', 0, 6460, -6366))).toBe('\u22126.4 kb')
  expect(markTag(mark('replacement', 0, 6000, 2000))).toBe('6 kb repl')
  expect(markTag(mark('replacement', 0, 176_000, 2800))).toBe('+2.8 kb')
  expect(markTag(mark('inversion', 0, 3000, 0))).toBe('3 kb inv')
  expect(
    markTag({ ...mark('duplication', 0, 0, 49_000), repeatedBp: 49_000 }),
  ).toBe('49 kb dup')
  expect(markTag({ ...mark('shift', 0, 5000, 0), shiftedBp: 2000 })).toBe(
    '2 kb shift',
  )
})

test('formTag names the largest change of a form and counts the rest', () => {
  expect(formTag({ states: [], marks: [] })).toBe('reference')
  expect(
    formTag({
      states: [NO_CALL_STATE, RARER_STATE],
      marks: [
        mark('deletion', 0, 2000, -2000),
        mark('deletion', 3000, 21_445, -18_445),
        mark('insertion', 30_000, 30_000, 300),
      ],
    }),
  ).toBe('18 kb del +4')
  expect(formTag({ states: [RARER_STATE], marks: [] })).toBe('rarer change')
  expect(formTag({ states: [NO_CALL_STATE], marks: [] })).toBe('not aligned')
})

test('formTag given the result leaves out the marks every form shares', () => {
  const shared = dep(1000, 3000, 0)
  const result = forms(
    [
      ...haps('R', 4, [shared]),
      ...haps('D', 4, [shared, dep(60_000, 70_000, 0)]),
    ],
    { minCarriers: 2 },
  )
  expect(result.sites.map(s => s.shared)).toEqual([true, false])
  expect(result.forms.map(f => formTag(f))).toEqual([
    '2 kb del',
    '10 kb del +1',
  ])
  expect(result.forms.map(f => formTag(f, result))).toEqual([
    'shared change',
    '10 kb del',
  ])
  const uncalled = { ...result.forms[0]!, states: [NO_CALL_STATE] }
  expect(formTag(uncalled, result)).toBe('not aligned')
})

test('a site most walks do not reach is an edge, left out of the key', () => {
  const result = forms(
    [
      ...Array.from({ length: 6 }, (_, i) => hap(`S${i}#1`, [], [[0, 90_000]])),
      ...haps('E', 3, [dep(95_000, 99_000, 0)]),
      ...haps('F', 2),
    ],
    { minCarriers: 2 },
  )
  expect(result.sites.map(s => [s.edge, s.covered])).toEqual([[true, 5]])
  expect(result.informative).toEqual([])
  expect(memberCounts(result)).toEqual([11])
})

test('forms a higher floor cannot tell apart merge', () => {
  const input = [
    ...haps('R', 5),
    ...haps('S', 5, [dep(1000, 1300, 0)]),
    ...haps('D', 5, [dep(5000, 25_000, 0)]),
  ]
  expect(memberCounts(forms(input, { floorBp: 50 }))).toEqual([5, 5, 5])
  expect(memberCounts(forms(input, { floorBp: 1000 }))).toEqual([10, 5])
})

test('the auto floor is about 10 px of the window, snapped to a preset', () => {
  expect(autoFloor(5000)).toBe(50)
  expect(autoFloor(24_000)).toBe(50)
  expect(autoFloor(35_000)).toBe(1000)
  expect(autoFloor(299_000)).toBe(1000)
  expect(autoFloor(500_000)).toBe(10_000)
  const input = [...haps('R', 2), ...haps('S', 2, [dep(1000, 1300, 0)])]
  expect(
    structuralForms(input, { window: { start: 0, end: 5000 } }).floorBp,
  ).toBe(50)
  expect(structuralForms(input, { window: WINDOW }).forms).toHaveLength(1)
})

test('excluded samples have no form, and none are excluded by default', () => {
  const input = [
    hap('CHM13#0', [dep(100, 1800, 0)]),
    hap('HG002#1'),
    hap('HG002#2'),
    hap('HG003#1'),
  ]
  expect(
    forms(input, { exclude: undefined, minCarriers: 1 })
      .forms.flatMap(f => f.members)
      .sort(),
  ).toEqual(['CHM13#0', 'HG002#1', 'HG002#2', 'HG003#1'])
  const result = forms(input, { exclude: ['CHM13'] })
  expect(result.forms.map(f => f.members)).toEqual([
    ['HG002#1', 'HG002#2', 'HG003#1'],
  ])
  expect(
    forms([hap('HG002#1'), hap('HG003#1')], { exclude: ['HG002'] }).forms[0]!
      .members,
  ).toEqual(['HG003#1'])
})

test('the result is plain data', () => {
  const result = forms([
    ...haps('R', 3),
    ...haps('D', 3, [dep(1000, 21_000, 0)]),
  ])
  expect(JSON.parse(JSON.stringify(result))).toEqual(result)
})

// GSTM1 at 462 haplotypes: 306 carry the exact 18,445 bp deletion, the rest
// R's structure, with SNP bubbles and indels under 50 bp throughout and two
// carriers of a 400 bp insertion too rare to tell apart
test('a GSTM1-like window is two forms', () => {
  const snps = (i: number) =>
    Array.from({ length: 6 }, (_, k) => {
      const at = 2000 + k * 5000 + (i % 7) * 13
      return dep(at, at + 1, 1)
    })
  const indel = (i: number) => (i % 5 === 0 ? [dep(30_100, 30_120, 0)] : [])
  const deletion = dep(7001, 25_447, 1)
  const input = Array.from({ length: 462 }, (_, i) => {
    const own = [...snps(i), ...indel(i)]
    if (i < 306) {
      own.push(deletion)
    }
    if (i === 400 || i === 401) {
      own.push(dep(33_000, 33_000, 400))
    }
    return hap(
      `HG${String(i).padStart(3, '0')}#1`,
      own.sort((a, b) => a.refStart - b.refStart),
      [[0, 35_000]],
    )
  })
  for (const floorBp of [50, 'auto'] as const) {
    const result = structuralForms(input, {
      window: { start: 0, end: 35_000 },
      floorBp,
    })
    expect(result.forms.map(f => [f.members.length, describeForm(f)])).toEqual([
      [306, '18 kb deletion'],
      [156, 'as the reference'],
    ])
    expect(result.minCarriers).toBe(5)
  }
  const at50 = structuralForms(input, {
    window: { start: 0, end: 35_000 },
    floorBp: 50,
  })
  expect(at50.rareCarriers).toEqual(['HG400#1', 'HG401#1'])
  expect(at50.forms[0]!.marks).toEqual([
    { ...mark('deletion', 7001, 25_447, -18_445), site: 0, state: '-18000' },
  ])
})
