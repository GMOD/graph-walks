import { TabixIndexedFile } from '@gmod/tabix'

import { RegionTooLargeError } from './util.ts'
import {
  ExpansionTooLargeError,
  formsCut,
  walkDepartures,
} from './walkDepartures.ts'
import { WalkReader } from './walkReader.ts'
import { walkCut } from './walkRows.ts'

import type { WalkDeparture } from './walkDeparture.ts'
import type { HaplotypeDepartures } from './walkDepartures.ts'
import type { Abridged, WalkNode, WalkRun } from './walkRows.ts'

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

const refStartOf = (nodes: Map<number, WalkNode>, id: number) => {
  const node = nodes.get(id)
  return node?.rank === 0 ? node.start : NaN
}

// A walk's run through the cut: `+3` steps forward through node 3, `-3`
// against it, from `hapStart` on its contig. R holds every node forward and
// once unless `against` or `refStarts` say otherwise
function walk(
  name: string,
  steps: string,
  more: {
    hapStart?: number
    abridged?: Abridged
    against?: number[]
    refStarts?: number[]
  } = {},
): WalkRun {
  const nodes = graphNodes()
  const parsed = steps
    .split(' ')
    .map(s => [+s.slice(1), s.startsWith('-')] as const)
  return {
    name,
    fragStart: 0,
    hapOffset: more.hapStart ?? 0,
    firstPiece: 0,
    pieceStarts: [],
    pieceChunks: [0],
    ids: Int32Array.from(parsed, ([id]) => id),
    rev: Uint8Array.from(parsed, ([, rev]) => (rev ? 1 : 0)),
    refStarts:
      more.refStarts === undefined
        ? Float64Array.from(parsed, ([id]) => refStartOf(nodes, id))
        : Float64Array.from(more.refStarts),
    against:
      more.against === undefined
        ? Uint8Array.from(parsed, ([, rev]) => (rev ? 1 : 0))
        : Uint8Array.from(more.against),
    ...(more.abridged ? { abridged: more.abridged } : {}),
  }
}

// A record of `name`'s departure: `ownBp` of its own from `hapStart` on its
// contig, between the leave and rejoin breakpoints, `+700` along R at 700
// and `-600` against it
function record(
  name: string,
  hapStart: number,
  ownBp: number,
  leave: string,
  rejoin: string,
  missingBp = 0,
): WalkDeparture {
  return {
    name,
    hapStart,
    ownBp,
    missingBp,
    leavePos: +leave.slice(1),
    leaveReversed: leave.startsWith('-'),
    rejoinPos: +rejoin.slice(1),
    rejoinReversed: rejoin.startsWith('-'),
    fromSteps: false,
  }
}

// The cut of these runs over the whole of R, the departures inside each run
// computed from its steps and `records` those between runs, read for forms
function departuresOf(runs: WalkRun[], records: WalkDeparture[] = []) {
  const nodes = graphNodes()
  const { fragments, departures } = walkCut(runs, nodes, R, records)
  return walkDepartures({ fragments, departures, nodes }, R)
}

const plain = {
  invertedBp: 0,
  repeatedBp: 0,
}

const detail = { refName: `${R}~0`, start: 0 }

test('a walk along R departs nowhere and covers what it reaches', () => {
  expect(departuresOf([walk('A#1#c', '+1 +2 +3')])).toEqual([
    { haplotype: 'A#1', covered: [[0, 700]], departures: [] },
  ])
})

test('a skipped stretch of R is a deletion with no ALT', () => {
  const [a] = departuresOf([walk('A#1#c', '+1 +3')])
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 600, altLength: 0, ...plain },
  ])
})

test('steps off R between abutting R steps insert at a point', () => {
  const [a] = departuresOf([walk('A#1#c', '+1 +50 +2')])
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 100, altLength: 300, ...plain },
  ])
})

test('steps off R in place of a stretch of R replace it', () => {
  const [a] = departuresOf([walk('A#1#c', '+1 +50 +3')])
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 600, altLength: 300, ...plain },
  ])
})

test('a stretch of R walked against it is ALT and inverted', () => {
  const [a] = departuresOf([walk('A#1#c', '+2 -3 +4')])
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
  const [a] = departuresOf([walk('A#1#c', '+1 +2 +3 +2 +3 +4')])
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
  const [a] = departuresOf([walk('A#1#c', '-4 -50 -2 -1')])
  expect(a!.departures).toEqual([
    { refStart: 600, refEnd: 700, altLength: 300, ...plain },
  ])
  expect(a!.covered).toEqual([[0, 1700]])
})

test('steps off R past either end of a walk are no departure', () => {
  expect(departuresOf([walk('A#1#c', '+51 +1 +2 +52')])).toEqual([
    { haplotype: 'A#1', covered: [[0, 600]], departures: [] },
  ])
  expect(departuresOf([walk('B#1#c', '+50 +52')])).toEqual([])
})

test('an abridged walk is its net change over its anchors on R', () => {
  const abridged = { pieces: [0], steps: 4000, bp: 1200, detail }
  const [a] = departuresOf([walk('A#1#c', '+2 +5', { abridged })])
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
  })
})

test("a haplotype's paths join into one list, the reference's are left out", () => {
  const result = departuresOf([
    walk(R, '+1 +2 +3 +4 +5 +6'),
    walk('A#1#c1', '+4 +6'),
    walk('A#1#c2', '+1 +50 +2'),
    walk('A#2#c1', '+1 +2'),
  ])
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

test('R that two paths of a haplotype cover, or a junction between them, has no call', () => {
  const [overlapping, abutting] = departuresOf([
    walk('A#1#c1', '+1 +2 +3 +4'),
    walk('A#1#c2', '+3 +4 +5'),
    walk('A#2#c1', '+1 +2'),
    walk('A#2#c2', '+3 +4'),
  ])
  expect(overlapping!.covered).toEqual([
    [0, 600],
    [1700, 1800],
  ])
  expect(abutting!.covered).toEqual([
    [0, 600],
    [600, 1700],
  ])
})

test('a departure with no step on the window is read from its record alone', () => {
  expect(departuresOf([], [record('A#1#c', 500, 0, '+100', '+1700')])).toEqual([
    {
      haplotype: 'A#1',
      covered: [[100, 1700]],
      departures: [{ refStart: 100, refEnd: 1700, altLength: 0, ...plain }],
    },
  ])
  expect(departuresOf([], [record('B#1#c', 500, 0, '-700', '-100')])).toEqual([
    {
      haplotype: 'B#1',
      covered: [[100, 700]],
      departures: [{ refStart: 100, refEnd: 700, altLength: 0, ...plain }],
    },
  ])
})

test("a departure across a contig's W lines counts the sequence the graph lacks as its own", () => {
  const [a] = departuresOf(
    [walk('A#1#c', '+1 +2'), walk('A#1#c', '+3 +4', { hapStart: 5600 })],
    [record('A#1#c', 600, 5000, '+600', '+600', 5000)],
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[0, 1700]],
    departures: [{ refStart: 600, refEnd: 600, altLength: 5000, ...plain }],
  })
})

test('an abridged fragment is net between its handles, and the departures either side stay exact', () => {
  const abridged = { pieces: [1], steps: 900, bp: 1500, detail }
  const [a] = departuresOf(
    [
      walk('A#1#c', '+1 +51'),
      walk('A#1#c', '+2 +4', { hapStart: 120, abridged }),
      walk('A#1#c', '+51 +5 +6', { hapStart: 1620 }),
    ],
    [
      record('A#1#c', 100, 20, '+100', '+100'),
      record('A#1#c', 1620, 20, '+1700', '+1700'),
    ],
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[0, 1900]],
    departures: [
      { refStart: 100, refEnd: 100, altLength: 20, ...plain },
      {
        refStart: 100,
        refEnd: 1700,
        altLength: 1500,
        ...plain,
        abridged: true,
      },
      { refStart: 1700, refEnd: 1700, altLength: 20, ...plain },
    ],
  })
})

test('a departure within an abridged stretch joins its net', () => {
  const abridged = { pieces: [1, 3], steps: 900, bp: 1200, detail }
  const [a] = departuresOf(
    [
      walk('A#1#c', '+2 +4', { hapStart: 100, abridged }),
      walk('A#1#c', '+3', { hapStart: 600 }),
    ],
    [
      record('A#1#c', 600, 0, '+600', '+600'),
      record('A#1#c', 700, 0, '+700', '+700'),
    ],
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[100, 1700]],
    departures: [
      {
        refStart: 100,
        refEnd: 1700,
        altLength: 1300,
        ...plain,
        abridged: true,
      },
    ],
  })
})

test('a walk along R through a node R holds reversed stays on R', () => {
  expect(
    departuresOf([walk('A#1#c', '-2 +3 +4', { against: [0, 0, 0] })]),
  ).toEqual([{ haplotype: 'A#1', covered: [[100, 1700]], departures: [] }])
  const [a] = departuresOf([walk('A#1#c', '+1 +50 -3', { against: [0, 0, 0] })])
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[0, 700]],
    departures: [{ refStart: 100, refEnd: 600, altLength: 300, ...plain }],
  })
})

test("a fragment's ends lie where their pass does, where R passes a node again", () => {
  const [a] = departuresOf(
    [walk('A#1#c', '+2 +4', { hapStart: 400, refStarts: [700, 1200] })],
    [record('A#1#c', 100, 300, '+100', '+700')],
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[100, 2200]],
    departures: [{ refStart: 100, refEnd: 700, altLength: 300, ...plain }],
  })
})

test("abridged pieces with another chunk's piece between them count its bp in their net", () => {
  const abridged = { pieces: [0, 2], steps: 3, bp: 1400, detail }
  const [a] = departuresOf(
    [
      walk('A#1#c', '+1 +4', { abridged }),
      walk('A#1#c', '+2 +3', { hapStart: 100 }),
      walk('A#1#c', '+5 +6', { hapStart: 2000 }),
    ],
    [record('A#1#c', 1700, 300, '+1700', '+1700')],
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[0, 1900]],
    departures: [
      {
        refStart: 0,
        refEnd: 1700,
        altLength: 2000,
        ...plain,
        abridged: true,
      },
    ],
  })
})

test("a departure leaving abridged pieces' last step keeps its own bp whole", () => {
  const abridged = { pieces: [0], steps: 2, bp: 600, detail }
  const [a] = departuresOf(
    [walk('A#1#c', '+1 +2', { abridged })],
    [record('A#1#c', 600, 5020, '+600', '+50000', 5000)],
  )
  expect(a!.departures).toEqual([
    { refStart: 0, refEnd: 600, altLength: 600, ...plain, abridged: true },
    { refStart: 600, refEnd: 50000, altLength: 5020, ...plain },
  ])
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

  test('with partial: false reads every abridged walk, or throws before reading any', async () => {
    const whole = byHaplotype(await formsCut(reader('unabridged.R'), region))
    expect(
      byHaplotype(
        await formsCut(reader('abridged.R'), region, { partial: false }),
      ),
    ).toEqual(whole)
    const error: unknown = await formsCut(reader('abridged.R'), region, {
      partial: false,
      stepBudget: 3,
    }).catch((e: unknown) => e)
    expect(error).toBeInstanceOf(ExpansionTooLargeError)
    expect(error).toBeInstanceOf(RegionTooLargeError)
    const { steps, limit, fitsBp } = error as ExpansionTooLargeError
    expect(limit).toBe(3)
    expect(steps).toBeGreaterThan(3)
    expect(fitsBp).toBe(Math.max(1, Math.floor((30 * 3) / steps)))
    expect(
      byHaplotype(
        await formsCut(reader('unabridged.R'), region, {
          partial: false,
          stepBudget: 0,
        }),
      ),
    ).toEqual(whole)
  })
})
