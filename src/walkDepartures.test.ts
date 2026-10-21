import { TabixIndexedFile } from '@gmod/tabix'

import { RegionTooLargeError } from './util.ts'
import { formsCut, walkDepartures } from './walkDepartures.ts'
import { WalkReader } from './walkReader.ts'

import type { GraphTables } from './util.ts'
import type { HaplotypeDepartures } from './walkDepartures.ts'
import type { WalkFragment, WalkNode } from './walkRows.ts'

const R = 'GRCh38#0#chr1'

// R as nodes 1-6 of 100, 500, 100, 1000, 100 and 100 bp from 0, and nodes 50
// and up off R
const refLengths = [100, 500, 100, 1000, 100, 100]
const offR: Record<number, number> = { 50: 300, 51: 20, 52: 2000 }

function graphNodes() {
  const nodes = new Map<number, WalkNode>()
  let start = 0
  refLengths.forEach((length, i) => {
    nodes.set(i + 1, {
      refName: R,
      start,
      end: start + length,
      rank: 0,
      onWindow: true,
    })
    start += length
  })
  for (const [id, length] of Object.entries(offR)) {
    nodes.set(+id, {
      refName: 'off',
      start: 0,
      end: length,
      rank: 1,
      onWindow: false,
    })
  }
  return nodes
}

// `+3` steps forward through node 3, `-3` against it
function fragment(
  name: string,
  steps: string,
  more: Partial<WalkFragment> = {},
): WalkFragment {
  const parsed = steps
    .split(' ')
    .map(s => [+s.slice(1), s.startsWith('-')] as const)
  return {
    name,
    fragStart: 0,
    hapStart: 0,
    hapEnd: 0,
    ids: Int32Array.from(parsed, ([id]) => id),
    rev: Uint8Array.from(parsed, ([, rev]) => (rev ? 1 : 0)),
    before: { kind: 'end' },
    after: { kind: 'end' },
    ...more,
  }
}

const departuresOf = (...fragments: WalkFragment[]) =>
  walkDepartures({ fragments, nodes: graphNodes() }, R)

const plain = {
  invertedBp: 0,
  repeatedBp: 0,
}

test('a walk along R departs nowhere and covers what it reaches', () => {
  expect(departuresOf(fragment('A#1#c', '+1 +2 +3'))).toEqual([
    { haplotype: 'A#1', covered: [[0, 700]], departures: [], open: false },
  ])
})

test('a skipped stretch of R is a deletion with no ALT', () => {
  const [a] = departuresOf(fragment('A#1#c', '+1 +3'))
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 600, altLength: 0, ...plain },
  ])
})

test('steps off R between abutting R steps insert at a point', () => {
  const [a] = departuresOf(fragment('A#1#c', '+1 +50 +2'))
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 100, altLength: 300, ...plain },
  ])
})

test('steps off R in place of a stretch of R replace it', () => {
  const [a] = departuresOf(fragment('A#1#c', '+1 +50 +3'))
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 600, altLength: 300, ...plain },
  ])
})

test('a stretch of R walked against it is ALT and inverted', () => {
  const [a] = departuresOf(fragment('A#1#c', '+2 -3 +4'))
  expect(a!.departures).toEqual([
    {
      refStart: 600,
      refEnd: 700,
      altLength: 100,
      invertedBp: 100,
      repeatedBp: 0,
    },
  ])
})

test('a step back along R is a tandem copy at the junction', () => {
  const [a] = departuresOf(fragment('A#1#c', '+1 +2 +3 +2 +3 +4'))
  expect(a!.departures).toEqual([
    {
      refStart: 700,
      refEnd: 700,
      altLength: 600,
      invertedBp: 0,
      repeatedBp: 600,
      copied: [100, 700],
    },
  ])
  expect(a!.covered).toEqual([[0, 1700]])
})

test('a walk reading R mostly reversed is read the other way round', () => {
  const [a] = departuresOf(fragment('A#1#c', '-4 -50 -2 -1'))
  expect(a!.departures).toEqual([
    { refStart: 600, refEnd: 700, altLength: 300, ...plain },
  ])
  expect(a!.covered).toEqual([[0, 1700]])
})

test('an abridged walk is its net change over its anchors on R', () => {
  const [a] = departuresOf(
    fragment('A#1#c', '+2 +5', {
      abridged: {
        pieces: [0],
        steps: 4000,
        bp: 1200,
        detail: { refName: `${R}~0`, start: 0 },
      },
    }),
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[100, 1800]],
    departures: [
      {
        refStart: 100,
        refEnd: 1800,
        altLength: 1200,
        ...plain,
        abridged: true,
      },
    ],
    open: false,
  })
})

test('steps off R past either end leave the walk open there', () => {
  const [a] = departuresOf(fragment('A#1#c', '+51 +1 +2 +52'))
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[0, 600]],
    departures: [],
    open: true,
  })
  const [b] = departuresOf(fragment('B#1#c', '+50 +52'))
  expect(b).toEqual({
    haplotype: 'B#1',
    covered: [],
    departures: [],
    open: true,
  })
})

test("a haplotype's fragments join into one list, the reference's are left out", () => {
  const result = departuresOf(
    fragment(R, '+1 +2 +3 +4 +5 +6'),
    fragment('A#1#c1', '+4 +6'),
    fragment('A#1#c2', '+1 +50 +2'),
    fragment('A#2#c1', '+1 +2'),
  )
  expect(
    result.map(h => [h.haplotype, h.covered, h.departures.length]),
  ).toEqual([
    [
      'A#1',
      [
        [0, 600],
        [700, 1900],
      ],
      2,
    ],
    ['A#2', [[0, 600]], 0],
  ])
  expect(result[0]!.departures.map(d => [d.refStart, d.refEnd])).toEqual([
    [100, 100],
    [1700, 1800],
  ])
})

test('tables give the same departures as fragments', () => {
  const names = ['1', '2', '3', '50']
  const tables: GraphTables = {
    nodes: {
      names,
      lengths: Int32Array.from([100, 500, 100, 300]),
      refs: Int32Array.from([0, 0, 0, -1]),
      starts: Float64Array.from([0, 100, 600, 0]),
      ranks: Int32Array.from([0, 0, 0, 1]),
      refNames: [R],
    },
    links: {
      from: new Int32Array(),
      to: new Int32Array(),
      strands: new Uint8Array(),
    },
    walks: {
      names: ['A#1#c', 'B#1#c'],
      starts: Float64Array.from([0, 0]),
      ends: Float64Array.from([0, 0]),
      offsets: Int32Array.from([0, 3, 5]),
      steps: Int32Array.from([0, 3, 2, 0, 2]),
      reversed: Uint8Array.from([0, 0, 0, 0, 0]),
      abridgedSteps: Float64Array.from([NaN, 900]),
      abridgedBp: Float64Array.from([NaN, 150]),
    },
  }
  expect(walkDepartures(tables, R).map(h => h.departures)).toEqual([
    [{ refStart: 100, refEnd: 600, altLength: 300, ...plain }],
    [
      {
        refStart: 0,
        refEnd: 700,
        altLength: 150,
        ...plain,
        abridged: true,
      },
    ],
  ])
})

test('an abridged fragment between two others is net between their R steps', () => {
  const abridged = {
    pieces: [1],
    steps: 900,
    bp: 1500,
    detail: { refName: `${R}~0`, start: 0 },
  }
  const [a] = departuresOf(
    fragment('A#1#c', '+1 +51', { hapStart: 0, hapEnd: 120 }),
    fragment('A#1#c', '+2 +4', { hapStart: 120, hapEnd: 1620, abridged }),
    fragment('A#1#c', '+51 +5 +6', { hapStart: 1620, hapEnd: 1840 }),
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[0, 1900]],
    departures: [
      {
        refStart: 100,
        refEnd: 1700,
        altLength: 20 + 1500 + 20,
        ...plain,
        abridged: true,
      },
    ],
    open: false,
  })
})

describe('a forms cut', () => {
  const DATA = new URL('../test_data/', import.meta.url).pathname
  const reader = (name: string) => {
    const path = `${DATA}${name}.walks.bed.gz`
    return new WalkReader(
      new TabixIndexedFile({ path, tbiPath: `${path}.tbi` }),
    )
  }
  const region = { refName: 'chr', start: 0, end: 30 }
  const byHaplotype = (hs: HaplotypeDepartures[]) =>
    Object.fromEntries(hs.map(h => [h.haplotype, h]))
  const net = (h: HaplotypeDepartures) =>
    h.departures.reduce(
      (sum, d) => sum + d.altLength - (d.refEnd - d.refStart),
      0,
    )

  test('reads abridged walks in full, as the unabridged file gives them', async () => {
    const whole = byHaplotype(await formsCut(reader('unabridged.R'), region))
    expect(byHaplotype(await formsCut(reader('abridged.R'), region))).toEqual(
      whole,
    )
    expect(Object.keys(whole).sort()).toEqual(['D#1', 'H#1', 'I#1', 'Q#0'])
  })

  test('keeps the net of each abridged walk past its step budget', async () => {
    const whole = await formsCut(reader('unabridged.R'), region)
    const abridged = byHaplotype(
      await formsCut(reader('abridged.R'), region, { stepBudget: 0 }),
    )
    for (const h of whole) {
      expect(abridged[h.haplotype]!.departures.every(d => d.abridged)).toBe(
        true,
      )
      expect(net(abridged[h.haplotype]!)).toBe(net(h))
    }
  })

  test('throws RegionTooLargeError past its byte budget', async () => {
    await expect(
      formsCut(reader('abridged.R'), region, { byteBudget: 1 }),
    ).rejects.toBeInstanceOf(RegionTooLargeError)
    expect(
      await formsCut(reader('abridged.R'), region, { byteBudget: 1e9 }),
    ).toEqual(await formsCut(reader('abridged.R'), region))
  })
})
