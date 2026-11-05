import { TabixIndexedFile } from '@gmod/tabix'

import { RegionTooLargeError } from './util.ts'
import {
  ExpansionTooLargeError,
  formsCut,
  walkDepartures,
} from './walkDepartures.ts'
import { DepartureReadError, WalkReader } from './walkReader.ts'
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
      passes: 1,
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
      passes: 0,
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
    leavePasses: 1,
    rejoinPos: +rejoin.slice(1),
    rejoinReversed: rejoin.startsWith('-'),
    rejoinPasses: 1,
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
  const abridged = { steps: 4000, bp: 1200, detail }
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
  const abridged = { steps: 900, bp: 1500, detail }
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

const netOf = (h: HaplotypeDepartures) =>
  h.departures.reduce(
    (sum, d) => sum + d.altLength - (d.refEnd - d.refStart),
    0,
  )

// What a path's abridged reading and its expanded one must agree on: each
// haplotype's net change and the R it covers
const netAndCover = (hs: HaplotypeDepartures[]) =>
  hs.map(h => ({ haplotype: h.haplotype, net: netOf(h), covered: h.covered }))

// A contig's W lines: an abridged piece over 2 and 3 running 100 bp of its
// own after them to its line's end, 50 bp no line covers, and an abridged
// piece of the next line running 100 bp of its own before 4
test('a departure across W lines leaving and rejoining inside abridged pieces is exact, and each piece a net without its own bp', () => {
  const [a] = departuresOf(
    [
      walk('A#1#c', '+2 +3', {
        hapStart: 100,
        abridged: { steps: 3, bp: 700, detail },
      }),
      walk('A#1#c', '+4 +4', {
        hapStart: 850,
        abridged: { steps: 2, bp: 1100, detail },
      }),
      walk('A#1#c', '+5 +6', { hapStart: 1950 }),
    ],
    [record('A#1#c', 700, 250, '+700', '+700', 50)],
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[100, 1900]],
    departures: [
      { refStart: 100, refEnd: 700, altLength: 600, ...plain, abridged: true },
      { refStart: 700, refEnd: 700, altLength: 250, ...plain },
      {
        refStart: 700,
        refEnd: 1700,
        altLength: 1000,
        ...plain,
        abridged: true,
      },
    ],
  })
})

test('a departure leaving inside an abridged piece keeps its own bp, and the piece ends at its leave', () => {
  const abridged = { steps: 900, bp: 1500, detail }
  const [a] = departuresOf(
    [
      walk('A#1#c', '+2 +4', { hapStart: 120, abridged }),
      walk('A#1#c', '+51 +5 +6', { hapStart: 1620 }),
    ],
    [record('A#1#c', 1600, 40, '+1700', '+1700')],
  )
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 1700, altLength: 1480, ...plain, abridged: true },
    { refStart: 1700, refEnd: 1700, altLength: 40, ...plain },
  ])
})

// A path 1 2 then 300 bp off R in node 50 and 2,000 in node 52, rejoining at
// 5: its first piece, 1 2 50, abridged, leaves R inside it, and 52 lies in a
// piece the cut does not read
test("a departure leaving inside an abridged piece's trailing own steps nets as the path read in full", () => {
  const unexpanded = departuresOf(
    [
      walk('A#1#c', '+1 +2', { abridged: { steps: 3, bp: 900, detail } }),
      walk('A#1#c', '+5 +6', { hapStart: 2900 }),
    ],
    [record('A#1#c', 600, 2300, '+600', '+1700')],
  )
  const expanded = departuresOf([walk('A#1#c', '+1 +2 +50 +52 +5 +6')])
  expect(netAndCover(expanded)).toEqual([
    { haplotype: 'A#1', net: 1200, covered: [[0, 1900]] },
  ])
  expect(netAndCover(unexpanded)).toEqual(netAndCover(expanded))
})

// The same path leaving R for good, rejoining it reversed past the window:
// its own bp lie on no stretch of R read here
test('a departure leaving inside an abridged piece and rejoining reversed nets as the path read in full', () => {
  const unexpanded = departuresOf(
    [walk('A#1#c', '+1 +2', { abridged: { steps: 3, bp: 900, detail } })],
    [record('A#1#c', 600, 2300, '+600', '-1700')],
  )
  const expanded = departuresOf(
    [walk('A#1#c', '+1 +2 +50 +52')],
    [record('A#1#c', 600, 2300, '+600', '-1700')],
  )
  expect(netAndCover(unexpanded)).toEqual(netAndCover(expanded))
})

// A path 1 2, 2,000 bp in node 52 and 300 in node 50, then 3 4 and 5 6: the
// piece 50 3 4 is abridged and the departure from 2 rejoins inside it
test("a departure rejoining inside an abridged piece's leading own steps nets as the path read in full", () => {
  const unexpanded = departuresOf(
    [
      walk('A#1#c', '+1 +2'),
      walk('A#1#c', '+3 +4', {
        hapStart: 2600,
        abridged: { steps: 3, bp: 1400, detail },
      }),
      walk('A#1#c', '+5 +6', { hapStart: 4000 }),
    ],
    [record('A#1#c', 600, 2300, '+600', '+600')],
  )
  const expanded = departuresOf([walk('A#1#c', '+1 +2 +52 +50 +3 +4 +5 +6')])
  expect(netAndCover(expanded)).toEqual([
    { haplotype: 'A#1', net: 2300, covered: [[0, 1900]] },
  ])
  expect(netAndCover(unexpanded)).toEqual(netAndCover(expanded))
})

// A piece 51 2 3 50, 20 bp of its own before R and 300 after, with a
// departure rejoining at its start's R and one leaving after 3, each running
// 2,000 bp more in node 52 outside it
test('departures at both ends of one abridged piece net as the path read in full', () => {
  const unexpanded = departuresOf(
    [
      walk('A#1#c', '+1'),
      walk('A#1#c', '+2 +3', {
        hapStart: 2100,
        abridged: { steps: 4, bp: 920, detail },
      }),
      walk('A#1#c', '+4 +5', { hapStart: 5020 }),
    ],
    [
      record('A#1#c', 100, 2020, '+100', '+100'),
      record('A#1#c', 2720, 2300, '+700', '+700'),
    ],
  )
  const expanded = departuresOf([
    walk('A#1#c', '+1 +52 +51 +2 +3 +50 +52 +4 +5'),
  ])
  expect(netAndCover(expanded)).toEqual([
    { haplotype: 'A#1', net: 4320, covered: [[0, 1800]] },
  ])
  expect(netAndCover(unexpanded)).toEqual(netAndCover(expanded))
})

test("a contig's second visit to the window gives the R both reach no call", () => {
  const [a] = departuresOf([
    walk('A#1#c', '+1 +2 +3'),
    walk('A#1#c', '+2 +3 +4', { hapStart: 5000 }),
  ])
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [
      [0, 100],
      [700, 1700],
    ],
    departures: [],
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

// An abridged row stands for one piece, from its hapOffset to hapOffset +
// bp, with its first and last steps on R as two steps, so pieces 0 and 2 abridged with piece 1 between them in another chunk
// are two nets, and piece 1 reads as itself between them
test('each abridged piece is a net over its own stretch, and a piece between two reads as itself', () => {
  const [a] = departuresOf([
    walk('A#1#c', '+1 +1', { abridged: { steps: 1, bp: 100, detail } }),
    walk('A#1#c', '+2 +3', { hapStart: 100 }),
    walk('A#1#c', '+4 +4', {
      hapStart: 700,
      abridged: { steps: 4, bp: 1300, detail },
    }),
    walk('A#1#c', '+5 +6', { hapStart: 2000 }),
  ])
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[0, 1900]],
    departures: [
      { refStart: 0, refEnd: 100, altLength: 100, ...plain, abridged: true },
      {
        refStart: 700,
        refEnd: 1700,
        altLength: 1300,
        ...plain,
        abridged: true,
      },
    ],
  })
})

test("a departure leaving abridged pieces' last step keeps its own bp whole", () => {
  const abridged = { steps: 2, bp: 600, detail }
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
      expect(abridged[h.haplotype]!.departures.some(d => d.abridged)).toBe(true)
      expect(net(abridged[h.haplotype]!)).toBe(net(h))
      expect(abridged[h.haplotype]!.covered).toEqual(h.covered)
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

  test("keeps an abridged walk's net exact when a piece of it lies in another chunk", async () => {
    // X#1's pieces 0 and 2 are abridged in chunk 0, and piece 1 between them
    // walks R in chunk 10
    const window = { refName: 'c', start: 0, end: 36 }
    const [abridged, expanded] = await Promise.all(
      [false, true].map(async expand =>
        byHaplotype(
          await formsCut(reader('interleaved.R'), window, { expand }),
        ),
      ),
    )
    const x = abridged!['X#1']!
    expect(net(expanded!['X#1']!)).toBe(84)
    expect(net(x)).toBe(84)
    expect(x.covered).toEqual(expanded!['X#1']!.covered)
  })

  // X#1 and Z#1 leave R inside an abridged piece whose last steps are their
  // own, V#1 and T#1 rejoin across W lines inside one whose first steps are
  // their own, Q#1 through a W line no piece holds, P#1 through a step on
  // another reference sequence, and S#1 rejoins reversed into a piece
  // against R
  test('reads each haplotype to the net and coverage its walks give in full when a departure leaves or rejoins inside an abridged piece, or a piece runs against R', async () => {
    const window = { refName: 'c', start: 0, end: 36 }
    const [abridged, expanded] = await Promise.all(
      [false, true].map(async expand =>
        (await formsCut(reader('ownends.R'), window, { expand }))
          .map(h => ({
            haplotype: h.haplotype,
            net: net(h),
            covered: h.covered,
            abridged: h.departures.some(d => d.abridged),
          }))
          .sort((a, b) => (a.haplotype < b.haplotype ? -1 : 1)),
      ),
    )
    expect(abridged!.every(h => h.abridged)).toBe(true)
    expect(abridged!.map(h => ({ ...h, abridged: false }))).toEqual(
      expanded!.map(h => ({ ...h, abridged: false })),
    )
    expect(expanded!.map(h => [h.haplotype, h.net])).toEqual([
      ['P#1', 70_110],
      ['Q#1', 70_028],
      ['S#1', 0],
      ['T#1', 70_008],
      ['V#1', 28],
      ['X#1', 18],
      ['Y#1', 40],
      ['Z#1', 10],
    ])
  })

  test('throws RegionTooLargeError past its byte budget', async () => {
    await expect(
      formsCut(reader('abridged.R'), region, { byteBudget: 1 }),
    ).rejects.toBeInstanceOf(RegionTooLargeError)
    expect(
      await formsCut(reader('abridged.R'), region, { byteBudget: 1e9 }),
    ).toEqual(await formsCut(reader('abridged.R'), region))
  })

  test('throws DepartureReadError past its departure byte budget', async () => {
    const window = { refName: 'chr1', start: 50, end: 55 }
    await expect(
      formsCut(reader('span.GRCh38'), window, { departureByteBudget: 1 }),
    ).rejects.toBeInstanceOf(DepartureReadError)
    expect(
      await formsCut(reader('span.GRCh38'), window, {
        departureByteBudget: 1e9,
      }),
    ).toEqual(await formsCut(reader('span.GRCh38'), window))
  })
})
