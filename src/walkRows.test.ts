import { graphTablesGFA } from './testing.ts'
import {
  WalkGraph,
  byteBudgetError,
  joinPieces,
  stepBudgetError,
  walkCut,
  walkHeader,
  walkNameFilter,
  walksKept,
} from './walkRows.ts'

import type { WalkDeparture } from './walkDeparture.ts'
import type { WalkNode, WalkRow } from './walkRows.ts'

// A decoded walk row: a piece of `name` from `hapOffset`, through `steps`,
// each placed on R where `nodes` holds its node on `refName`
function walkRow(
  name: string,
  piece: number,
  hapOffset: number,
  steps: [number, boolean][],
  more: Partial<WalkRow> = {},
  nodes?: Map<number, WalkNode>,
  refName = 'GRCh38#0#chr1',
): WalkRow {
  const onR = (id: number) => {
    const node = nodes?.get(id)
    return node?.rank === 0 && node.refName === refName ? node.start : NaN
  }
  return {
    chunkStart: 0,
    name,
    fragStart: 0,
    hapOffset,
    piece,
    n: steps.length,
    steps: {
      ids: Int32Array.from(steps, ([id]) => id),
      rev: Uint8Array.from(steps, ([, reversed]) => (reversed ? 1 : 0)),
      refStarts: Float64Array.from(steps, ([id]) => onR(id)),
      against: Uint8Array.from(steps, ([, reversed]) => (reversed ? 1 : 0)),
    },
    ...more,
  }
}

const forward = (...ids: number[]) =>
  ids.map(id => [id, false] as [number, boolean])

test('continuation rows join by piece index, and a missing piece splits the walk', () => {
  const rows = [
    walkRow('HG002#1#chr1', 3, 900, forward(30)),
    walkRow('HG002#1#chr1', 1, 300, forward(12, 13)),
    walkRow('HG002#1#chr1', 0, 100, forward(10, 11)),
  ]
  const runs = joinPieces(rows)
  expect(runs.map(r => [r.hapOffset, [...r.ids], r.pieceStarts])).toEqual([
    [100, [10, 11, 12, 13], [2]],
    [900, [30], []],
  ])
})

test("pieces join within their W line: a contig's second W line, numbering its pieces from 0 too, joins none of the first's", () => {
  const rows = [
    walkRow('H#1#h', 1, 20, forward(3, 4), { fragStart: 100 }),
    walkRow('H#1#h', 0, 0, forward(1, 2)),
    walkRow('H#1#h', 0, 0, forward(1, 2), { fragStart: 100 }),
    walkRow('H#1#h', 1, 20, [
      [3, true],
      [4, false],
    ]),
  ]
  expect(
    joinPieces(rows)
      .map(r => [r.fragStart, [...r.ids], [...r.rev]])
      .sort(([a], [b]) => Number(a) - Number(b)),
  ).toEqual([
    [0, [1, 2, 3, 4], [0, 0, 1, 0]],
    [100, [1, 2, 3, 4], [0, 0, 0, 0]],
  ])
})

test('the haplotype filter keeps the named walks and the reference', () => {
  expect(walkNameFilter(undefined, 'GRCh38#0')).toBeUndefined()
  expect(walkNameFilter([], 'GRCh38#0')).toBeUndefined()
  const keep = walkNameFilter(['HG002', 'HG00097#2'], 'GRCh38#0')!
  expect(
    [
      'GRCh38#0#chr1',
      'HG002#1#chr1',
      'HG002#2#chr1',
      'HG00097#1#CM1',
      'HG00097#2#CM2',
      'HG0029#1#chr1',
    ].filter(keep),
  ).toEqual(['GRCh38#0#chr1', 'HG002#1#chr1', 'HG002#2#chr1', 'HG00097#2#CM2'])
})

test("a cut keeps its reference walk by the reference's whole name when it has no haplotype", () => {
  const keep = walksKept(['HG002'], 'chm13#chr1')!
  expect(['chm13#chr1', 'HG002#1#chr1', 'HG003#1#chr1'].filter(keep)).toEqual([
    'chm13#chr1',
    'HG002#1#chr1',
  ])
  expect(walksKept(['HG002'], 'GRCh38#0#chr1')!('GRCh38#0#chr1')).toBe(true)
})

describe('the step budget', () => {
  const chunk = 100
  const rows = (perChunk: [number, number][]): WalkRow[] =>
    perChunk.map(([chunkStart, n]) => ({
      name: 'HG002#1#chr1',
      fragStart: 0,
      hapOffset: 0,
      piece: 0,
      n,
      interior: n,
      steps: { ids: new Int32Array(n), rev: new Uint8Array(n) },
      chunkStart,
    }))

  test('passes a window under it', () => {
    expect(
      stepBudgetError(
        rows([
          [0, 5],
          [100, 5],
        ]),
        10,
        1,
        { start: 150, end: 190 },
        chunk,
        false,
      ),
    ).toBeUndefined()
  })

  test('names the span from the window start that fits', () => {
    const error = stepBudgetError(
      rows([
        [100, 2],
        [200, 3],
        [300, 3],
        [400, 50],
      ]),
      10,
      1,
      { start: 250, end: 480 },
      chunk,
      false,
    )!
    expect(error.name).toBe('NodeLimitError')
    expect(error.regionTooLarge).toBe(true)
    expect(error.fitsBp).toBe(150)
    expect(error.zoomHelps).toBe(true)
    expect(error.bytes).toBeUndefined()
    expect(error.message).toBe('Zoom in to about 150bp to see the graph')
  })

  test('offers a cut for every haplotype fewer as well as a zoom', () => {
    // four walks, the reference's among them, with these steps each per chunk
    const each: [number, number][] = [
      [100, 1],
      [200, 1],
      [300, 2],
      [400, 1],
    ]
    const four = [
      'GRCh38#0#chr1',
      'HG002#1#chr1',
      'HG002#2#chr1',
      'NA19240#1#chr1',
    ].flatMap(name => rows(each).map(row => ({ ...row, name })))
    const window = { start: 250, end: 480 }
    expect(stepBudgetError(four, 10, 1, window, chunk, false)!.message).toBe(
      'Zoom in to about 50bp, or choose fewer haplotypes, to see the graph',
    )
    expect(stepBudgetError(four, 10, 1, window, chunk, true)!.message).toBe(
      'Zoom in to about 50bp to see the graph',
    )
  })

  test('asks for fewer haplotypes where no zoom fits', () => {
    const error = stepBudgetError(
      rows([
        [100, 50],
        [200, 3],
      ]),
      5,
      2,
      { start: 250, end: 280 },
      chunk,
      true,
    )!
    expect(error.message).toBe(
      'Too many haplotype steps here to draw these haplotypes (50 steps read, a step along the reference counting a quarter, over 2 × the step budget of 5); choose fewer haplotypes',
    )
  })

  test('counts a step along the reference as a quarter of one off it', () => {
    const along = rows([[100, 40]]).map(row => ({ ...row, interior: 0 }))
    expect(
      stepBudgetError(along, 5, 2, { start: 150, end: 180 }, chunk, true),
    ).toBe(undefined)
    expect(
      stepBudgetError(
        rows([[100, 40]]),
        5,
        2,
        { start: 150, end: 180 },
        chunk,
        true,
      ),
    ).toBeDefined()
  })
})

test("a walk file's header states its format, tile, references, haplotypes and holes", () => {
  expect(walkHeader(['#walks\ttile:i:65536'])).toEqual({
    tile: 65_536,
    references: [],
    haplotypes: [],
    holes: [],
  })
  expect(walkHeader(['#walks']).tile).toBeUndefined()
  expect(
    walkHeader([
      '#walks\tformat:i:7\ttile:i:1000\tstay:i:311',
      '#reference\tGRCh38',
      '#reference\tCHM13',
      '#haplotype\tHG00097#1',
      '#haplotype\tHG002#1',
      '#hole\tGRCh38#0#chr1\t12\t20',
    ]),
  ).toEqual({
    format: 7,
    tile: 1000,
    references: ['GRCh38', 'CHM13'],
    haplotypes: ['HG00097#1', 'HG002#1'],
    holes: [{ refName: 'GRCh38#0#chr1', start: 12, end: 20 }],
  })
})

describe('the byte budget', () => {
  const tile = 100
  // bytes a read from the window's first base to `end` fetches, 10 per tile
  // and 80 for the tile at 300
  const bytesTo = (end: number) => {
    let bytes = 0
    for (let cs = 0; cs < end; cs += tile) {
      bytes += cs === 300 ? 80 : 10
    }
    return bytes
  }

  test('passes a window under it', () => {
    expect(
      byteBudgetError(bytesTo, 40, { start: 150, end: 250 }, tile),
    ).toBeUndefined()
  })

  test('names the span from the window start that fits', () => {
    const error = byteBudgetError(bytesTo, 50, { start: 150, end: 450 }, tile)!
    expect(error.name).toBe('NodeLimitError')
    expect(error.regionTooLarge).toBe(true)
    expect(error.fitsBp).toBe(150)
    expect(error.zoomHelps).toBe(true)
    expect(error.bytes).toBe(120)
    expect(error.message).toBe('Zoom in to about 150bp to see the graph')
  })

  test('says so where no zoom fits', () => {
    const error = byteBudgetError(
      end => bytesTo(end) * 100_000,
      5_000_000,
      { start: 310, end: 320 },
      tile,
    )!
    expect(error.message).toMatch(
      /^Too much graph here to fetch, even zoomed in \(11\.0 MB against a byte budget of 5\.0 MB\)\. A satellite or rDNA array/,
    )
    expect(error.zoomHelps).toBe(false)
    expect(error.bytes).toBe(11_000_000)
  })
})

// Reference 1-2-3-4, 100 bp each from 0; the window is 150-250, so only 2 and
// 3 are on it. Allele 5 replaces 2 and allele 6 sits before 1.
function graph() {
  const g = new WalkGraph({
    refName: 'GRCh38#0#chr1',
    start: 150,
    end: 250,
    context: 0,
  })
  const node = (id: number, rank: number, refName: string, start: number) => {
    g.placeNode(id, {
      rank,
      refName,
      start,
      end: start + 100,
      passes: rank === 0 ? 1 : 0,
    })
  }
  node(1, 0, 'GRCh38#0#chr1', 0)
  node(2, 0, 'GRCh38#0#chr1', 100)
  node(3, 0, 'GRCh38#0#chr1', 200)
  node(4, 0, 'GRCh38#0#chr1', 300)
  node(5, 1, 'HG002#1#chr1', 1000)
  node(6, 1, 'HG002#1#chr1', 2000)
  node(2, 0, 'GRCh38#0#chr1', 100)
  for (const [a, b] of [
    [1, 2],
    [2, 3],
    [3, 4],
    [1, 5],
    [2, 5],
    [5, 3],
    [6, 1],
  ]) {
    g.addLinkEnds(a! * 2, b! * 2)
  }
  return g
}

// a run of `name` through `steps`, placed on the graph's reference
const run = (
  g: WalkGraph,
  name: string,
  hapOffset: number,
  steps: number[],
  more: Partial<WalkRow> = {},
) =>
  walkRow(
    name,
    0,
    hapOffset,
    steps.map(id => [id, false] as [number, boolean]),
    more,
    g.nodes,
  )

const departure = (d: WalkDeparture) => [
  d.name,
  d.hapStart,
  d.ownBp,
  `${d.leavePos}${d.leaveReversed ? '-' : '+'}`,
  `${d.rejoinPos}${d.rejoinReversed ? '-' : '+'}`,
  d.fromSteps ? 'steps' : 'record',
  d.steps ? [...d.steps.ids] : undefined,
]

test('a node is placed once', () => {
  const g = graph()
  expect(g.nodes.size).toBe(6)
  expect(g.nodes.get(2)).toEqual({
    rank: 0,
    refName: 'GRCh38#0#chr1',
    start: 100,
    end: 200,
    passes: 1,
    onWindow: true,
  })
})

test("a node's reference position wins, whichever chunk is read first", () => {
  const g = graph()
  g.placeNode(7, {
    rank: 1,
    refName: 'HG002#1#chr1',
    start: 0,
    end: 100,
    passes: 0,
  })
  g.placeNode(7, {
    rank: 0,
    refName: 'GRCh38#0#chr1',
    start: 200,
    end: 300,
    passes: 1,
  })
  g.placeNode(7, {
    rank: 1,
    refName: 'HG002#1#chr1',
    start: 0,
    end: 100,
    passes: 0,
  })
  g.placeNode(8, {
    rank: 0,
    refName: 'GRCh38#0#chr1',
    start: 900,
    end: 1000,
    passes: 1,
  })
  g.placeNode(8, {
    rank: 0,
    refName: 'GRCh38#0#chr1',
    start: 150,
    end: 250,
    passes: 1,
  })
  expect([g.nodes.get(7), g.nodes.get(8)]).toEqual([
    {
      rank: 0,
      refName: 'GRCh38#0#chr1',
      start: 200,
      end: 300,
      passes: 1,
      onWindow: true,
    },
    {
      rank: 0,
      refName: 'GRCh38#0#chr1',
      start: 150,
      end: 250,
      passes: 1,
      onWindow: true,
    },
  ])
})

test("a node's pass over the window wins over the reference's other passes, whichever is read first", () => {
  const g = graph()
  const at = (start: number) => ({
    rank: 0,
    refName: 'GRCh38#0#chr1',
    start,
    end: start + 100,
    passes: 3,
  })
  g.placeNode(9, at(0))
  g.placeNode(9, at(200))
  g.placeNode(9, at(160))
  g.placeNode(10, at(160))
  g.placeNode(10, at(0))
  expect([g.nodes.get(9), g.nodes.get(10)]).toEqual([
    { ...at(160), onWindow: true },
    { ...at(160), onWindow: true },
  ])
})

test.each([false, true])(
  'off the window, a node takes the reference row whose path sorts first, then its earliest start, whichever is read first: reversed %s',
  reversed => {
    const g = graph()
    const at = (refName: string, start: number) => ({
      rank: 0,
      refName,
      start,
      end: start + 100,
      passes: 2,
    })
    const placements = [
      at('GRCh38#0#chr1', 900),
      at('GRCh38#0#chr0', 1500),
      at('GRCh38#0#chr0', 1200),
    ]
    for (const node of reversed ? placements.reverse() : placements) {
      g.placeNode(11, node)
    }
    expect(g.nodes.get(11)).toEqual({
      ...at('GRCh38#0#chr0', 1200),
      onWindow: false,
    })
  },
)

test('a cut is written as the GFA its tables stand for, with the departures over it', () => {
  const g = graph()
  const runs = joinPieces([
    run(g, 'GRCh38#0#chr1', 0, [1, 2, 3, 4]),
    // brings allele 5 into the cut between two steps on the window
    run(g, 'HG002#1#chr1', 0, [2, 5, 3]),
    // enters at 3 from 1 through allele 5, in place of 2
    run(g, 'HG002#2#chr1', 1000, [6, 1, 5, 3, 4]),
  ])
  const { kept, fragments, departures } = walkCut(
    runs,
    g.nodes,
    'GRCh38#0#chr1',
    [],
  )
  expect([...kept].sort()).toEqual([2, 3, 5])
  expect(
    fragments.map(f => [
      f.name,
      f.hapStart,
      f.hapEnd,
      [...f.ids],
      f.trimmedBefore,
      f.trimmedAfter,
    ]),
  ).toEqual([
    ['GRCh38#0#chr1', 100, 300, [2, 3], 100, 100],
    ['HG002#1#chr1', 0, 300, [2, 5, 3], 0, 0],
    ['HG002#2#chr1', 1300, 1400, [3], 300, 100],
  ])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 100, 100, '200+', '200+', 'steps', [5]],
    ['HG002#2#chr1', 1200, 100, '100+', '200+', 'steps', [5]],
  ])
  const tables = g.tables(kept, fragments, departures)
  expect(graphTablesGFA(tables).split('\n')).toEqual([
    'H\tVN:Z:1.1',
    'S\t2\t*\tLN:i:100\tSN:Z:GRCh38#0#chr1\tSO:i:100\tSR:i:0',
    'S\t3\t*\tLN:i:100\tSN:Z:GRCh38#0#chr1\tSO:i:200\tSR:i:0',
    'S\t5\t*\tLN:i:100\tSN:Z:HG002#1#chr1\tSO:i:1000\tSR:i:1',
    'L\t2\t+\t3\t+\t0M',
    'L\t2\t+\t5\t+\t0M',
    'L\t5\t+\t3\t+\t0M',
    'W\tGRCh38\t0\tchr1\t100\t300\t>2>3',
    'W\tHG002\t1\tchr1\t0\t300\t>2>5>3',
    'W\tHG002\t2\tchr1\t1300\t1400\t>3',
  ])
  expect([...tables.walks.trimmedBefore!]).toEqual([100, 0, 300])
  expect([...tables.walks.trimmedAfter!]).toEqual([100, 0, 100])
  const d = tables.departures
  expect(d.names).toEqual(['HG002#1#chr1', 'HG002#2#chr1'])
  expect([...d.hapStarts]).toEqual([100, 1200])
  expect([...d.leavePos]).toEqual([200, 100])
  expect([...d.rejoinPos]).toEqual([200, 200])
  expect([...d.fromSteps]).toEqual([1, 1])
  expect([...d.offsets]).toEqual([0, 1, 2])
  expect([...d.steps].map(i => tables.nodes.names[i])).toEqual(['5', '5'])
})

test('a walk leaving the window is cut at its last step on it, and its departure says where it rejoins', () => {
  const g = graph()
  const runs = joinPieces([
    run(g, 'GRCh38#0#chr1', 0, [1, 2, 3, 4]),
    // leaves at 2 through allele 5 and rejoins at 4
    run(g, 'HG002#1#chr1', 0, [1, 2, 5, 4]),
    // from 2 straight to 4
    run(g, 'HG002#2#chr1', 0, [1, 2, 4]),
    // enters at 3 from 1 through allele 6
    run(g, 'HG002#3#chr1', 0, [5, 1, 6, 3]),
    // leaves at 2 through allele 5 and ends there
    run(g, 'HG002#4#chr1', 0, [2, 5]),
    // starts on allele 6 and enters at 3
    run(g, 'HG002#5#chr1', 0, [6, 3]),
  ])
  const { kept, fragments, departures } = walkCut(
    runs,
    g.nodes,
    'GRCh38#0#chr1',
    [],
  )
  expect([...kept].sort()).toEqual([2, 3, 5, 6])
  expect(
    fragments.map(f => [f.name, f.hapStart, f.hapEnd, [...f.ids]]),
  ).toEqual([
    ['GRCh38#0#chr1', 100, 300, [2, 3]],
    ['HG002#1#chr1', 100, 200, [2]],
    ['HG002#2#chr1', 100, 200, [2]],
    ['HG002#3#chr1', 300, 400, [3]],
    ['HG002#4#chr1', 0, 100, [2]],
    ['HG002#5#chr1', 100, 200, [3]],
  ])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 200, 100, '200+', '300+', 'steps', [5]],
    ['HG002#2#chr1', 200, 0, '200+', '300+', 'steps', undefined],
    ['HG002#3#chr1', 200, 100, '100+', '200+', 'steps', [6]],
  ])
})

test('a departure between the pieces of one run comes from the records, not the steps', () => {
  const g = graph()
  const rows = [
    run(g, 'HG002#1#chr1', 0, [1, 2, 5]),
    walkRow('HG002#1#chr1', 1, 300, forward(4), { chunkStart: 300 }, g.nodes),
  ]
  const record: WalkDeparture = {
    name: 'HG002#1#chr1',
    hapStart: 200,
    ownBp: 100,
    missingBp: 0,
    leavePos: 200,
    leaveReversed: false,
    leavePasses: 1,
    rejoinPos: 300,
    rejoinReversed: false,
    rejoinPasses: 1,
    fromSteps: false,
  }
  const { fragments, departures, unattached } = walkCut(
    joinPieces(rows),
    g.nodes,
    'GRCh38#0#chr1',
    [record],
  )
  expect(fragments.map(f => [f.name, [...f.ids], f.trimmedAfter])).toEqual([
    ['HG002#1#chr1', [2], 200],
  ])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 200, 100, '200+', '300+', 'record', [5]],
  ])
  expect(unattached).toEqual([])
})

test("a record's own steps come from its leave and rejoin pieces, and a record lacking either is asked for", () => {
  const g = graph()
  const record = (hapStart: number): WalkDeparture => ({
    name: 'HG002#1#chr1',
    hapStart,
    ownBp: 100,
    missingBp: 0,
    leavePos: 200,
    leaveReversed: false,
    leavePasses: 1,
    rejoinPos: 300,
    rejoinReversed: false,
    rejoinPasses: 1,
    fromSteps: false,
  })
  const rejoinOnly = walkRow(
    'HG002#1#chr1',
    1,
    300,
    forward(4),
    { chunkStart: 300 },
    g.nodes,
  )
  const lacking = walkCut(joinPieces([rejoinOnly]), g.nodes, 'GRCh38#0#chr1', [
    record(200),
  ])
  expect(lacking.departures.map(departure)).toEqual([
    ['HG002#1#chr1', 200, 100, '200+', '300+', 'record', undefined],
  ])
  expect(lacking.unattached).toEqual([{ departure: record(200) }])
  const leaveOnly = walkCut(
    joinPieces([run(g, 'HG002#1#chr1', 0, [1, 2])]),
    g.nodes,
    'GRCh38#0#chr1',
    [record(200)],
  )
  expect(leaveOnly.unattached).toEqual([{ departure: record(200) }])
  const over = walkCut(
    joinPieces([run(g, 'HG002#1#chr1', 0, [1, 2, 5]), rejoinOnly]),
    g.nodes,
    'GRCh38#0#chr1',
    [record(200)],
    { attachBp: 50 },
  )
  expect(over.departures.map(departure)).toEqual([
    ['HG002#1#chr1', 200, 100, '200+', '300+', 'record', undefined],
  ])
  expect(over.unattached).toEqual([])
})

test('a contained cut splits a walk where it steps onto the reference off the window', () => {
  const g = graph()
  g.placeNode(8, {
    rank: 1,
    refName: 'HG002#3#chr1',
    start: 0,
    end: 100,
    passes: 0,
  })
  const runs = joinPieces([
    run(g, 'GRCh38#0#chr1', 0, [1, 2, 3, 4]),
    // takes 2 and 3, goes back to 1 and takes 2 again, then 4
    run(g, 'HG002#1#chr1', 0, [2, 3, 1, 2, 4]),
    // leaves 3 through alleles 5 and 6 and comes back to 2
    run(g, 'HG002#2#chr1', 0, [2, 3, 5, 6, 2]),
    // leaves 3 through alleles 5 and 8 and does not come back
    run(g, 'HG002#3#chr1', 0, [3, 5, 8]),
  ])
  const { kept, fragments, departures } = walkCut(
    runs,
    g.nodes,
    'GRCh38#0#chr1',
    [],
    { contained: true, context: 50 },
  )
  expect([...kept].sort()).toEqual([2, 3, 5, 6])
  expect(
    fragments.map(f => [
      f.name,
      f.hapStart,
      f.hapEnd,
      [...f.ids],
      f.trimmedBefore,
      f.trimmedAfter,
    ]),
  ).toEqual([
    ['GRCh38#0#chr1', 100, 300, [2, 3], 100, 100],
    ['HG002#1#chr1', 0, 200, [2, 3], 0, 300],
    ['HG002#1#chr1', 300, 400, [2], 300, 100],
    ['HG002#2#chr1', 0, 500, [2, 3, 5, 6, 2], 0, 0],
    ['HG002#3#chr1', 0, 200, [3, 5], 0, 100],
  ])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 200, 0, '300+', '0+', 'steps', undefined],
    ['HG002#1#chr1', 400, 0, '200+', '300+', 'steps', undefined],
    ['HG002#2#chr1', 200, 200, '300+', '100+', 'steps', [5, 6]],
  ])
})

test('the reference walk is written first, whatever order the rows came in', () => {
  const g = graph()
  const runs = joinPieces([
    run(g, 'HG002#2#chr1', 0, [2, 3]),
    run(g, 'GRCh38#0#chr1', 0, [1, 2, 3, 4]),
    run(g, 'HG002#1#chr1', 0, [2, 5, 3]),
  ])
  const { kept, fragments, departures } = walkCut(
    runs,
    g.nodes,
    'GRCh38#0#chr1',
    [],
  )
  expect(g.tables(kept, fragments, departures).walks.names).toEqual([
    'GRCh38#0#chr1',
    'HG002#1#chr1',
    'HG002#2#chr1',
  ])
})

describe('the stretches between and past a walk’s fragments', () => {
  const record = (more: Partial<WalkDeparture>): WalkDeparture => ({
    name: 'HG002#1#chr1',
    hapStart: 200,
    ownBp: 300,
    missingBp: 200,
    leavePos: 200,
    leaveReversed: false,
    leavePasses: 1,
    rejoinPos: 200,
    rejoinReversed: false,
    rejoinPasses: 1,
    fromSteps: false,
    ...more,
  })
  const segments = (gaps: ReturnType<ReturnType<typeof walkCut>['gapsOf']>) =>
    gaps.map(g => [
      g.edge,
      g.start,
      g.end,
      g.segments.map(x =>
        x.kind === 'steps'
          ? [x.kind, x.bp, x.ids]
          : x.kind === 'own'
            ? [x.kind, x.bp, x.missingBp]
            : x.kind === 'unread'
              ? [x.kind, x.bp, x.before, x.after]
              : [x.kind, x.bp],
      ),
    ])
  // a W line that leaves at 2 for allele 5 and ends, 200 bp the graph lacks,
  // and the next W line rejoining at 3, as HG00097#1 does at the amylase locus
  const lineEnd = (g: WalkGraph) =>
    joinPieces([
      run(g, 'HG002#1#chr1', 100, [2, 5]),
      run(g, 'HG002#1#chr1', 500, [3], { fragStart: 500 }),
    ])

  test('a record’s attached steps, then its missing bp, fill the stretch between two W lines', () => {
    const g = graph()
    const { gapsOf } = walkCut(lineEnd(g), g.nodes, 'GRCh38#0#chr1', [
      record({}),
    ])
    expect(segments(gapsOf())).toEqual([
      [
        0,
        200,
        500,
        [
          ['steps', 100, [5]],
          ['own', 200, 200],
        ],
      ],
    ])
  })

  test('a record over attachBp leaves its steps read and kept out', () => {
    const g = graph()
    const { gapsOf } = walkCut(
      lineEnd(g),
      g.nodes,
      'GRCh38#0#chr1',
      [record({})],
      { attachBp: 50 },
    )
    expect(segments(gapsOf())).toEqual([
      [
        0,
        200,
        500,
        [
          ['read', 100],
          ['own', 200, 200],
        ],
      ],
    ])
  })

  test('pieces the cut did not read, with no record over them, are unread between the breakpoints either side', () => {
    const g = graph()
    const runs = joinPieces([
      run(g, 'HG002#1#chr1', 100, [2]),
      walkRow('HG002#1#chr1', 5, 1000, forward(3), {}, g.nodes),
    ])
    const { gapsOf } = walkCut(runs, g.nodes, 'GRCh38#0#chr1', [])
    expect(segments(gapsOf())).toEqual([
      [0, 200, 1000, [['unread', 800, 200, 200]]],
    ])
  })

  test('a row’s end reaches over a record leaving at it, and a record past the row stands alone', () => {
    const g = graph()
    const { gapsOf } = walkCut(
      joinPieces([run(g, 'HG002#1#chr1', 100, [2])]),
      g.nodes,
      'GRCh38#0#chr1',
      [
        record({ missingBp: 0, rejoinPos: 300 }),
        record({
          hapStart: 900,
          ownBp: 50,
          missingBp: 0,
          leavePos: 300,
          rejoinPos: 400,
        }),
      ],
    )
    expect(segments(gapsOf())).toEqual([
      [2, 200, 500, [['own', 300, 0]]],
      [3, 900, 950, [['own', 50, 0]]],
    ])
  })

  test('a walk with no fragment has a stretch for each of its records', () => {
    const g = graph()
    const { gapsOf } = walkCut([], g.nodes, 'GRCh38#0#chr1', [record({})])
    expect(segments(gapsOf())).toEqual([[3, 200, 500, [['own', 300, 200]]]])
  })

  test('the tables key each stretch to its fragment’s walk and its steps to node indices', () => {
    const g = graph()
    const { kept, fragments, departures, gapsOf } = walkCut(
      lineEnd(g),
      g.nodes,
      'GRCh38#0#chr1',
      [record({})],
    )
    const tables = g.tables(kept, fragments, departures, gapsOf())
    const t = tables.gaps!
    expect(tables.walks.names).toEqual(['HG002#1#chr1', 'HG002#1#chr1'])
    expect([...t.walks]).toEqual([0])
    expect([...t.segments.kinds]).toEqual([0, 3])
    expect([...t.segments.steps].map(i => tables.nodes.names[i])).toEqual(['5'])
  })

  test('a computed departure’s steps past the row’s last fragment are its steps, and the reference past them is read', () => {
    const g = graph()
    const { gapsOf } = walkCut(
      joinPieces([run(g, 'HG002#1#chr1', 100, [2, 5, 4])]),
      g.nodes,
      'GRCh38#0#chr1',
      [],
    )
    expect(segments(gapsOf())).toEqual([
      [
        2,
        200,
        400,
        [
          ['steps', 100, [5]],
          ['read', 100],
        ],
      ],
    ])
  })

  test('a tile fetched for a record’s own steps gives those steps and nothing read', () => {
    const g = graph()
    const { gapsOf } = walkCut(
      joinPieces([run(g, 'HG002#1#chr1', 100, [2, 5])]),
      g.nodes,
      'GRCh38#0#chr1',
      [record({ missingBp: 100 })],
      {
        attach: joinPieces([
          run(g, 'HG002#1#chr1', 400, [6, 3], { fragStart: 400 }),
        ]),
      },
    )
    expect(segments(gapsOf())).toEqual([
      [
        2,
        200,
        500,
        [
          ['steps', 100, [5]],
          ['own', 100, 100],
          ['steps', 100, [6]],
        ],
      ],
    ])
  })

  test('an abridged run between fragments is abridged, and its handles are where the walk was last and next seen', () => {
    const g = graph()
    for (const [id, start] of [
      [7, 400],
      [8, 500],
    ] as const) {
      g.placeNode(id, {
        rank: 0,
        refName: 'GRCh38#0#chr1',
        start,
        end: start + 100,
        passes: 1,
      })
    }
    const detail = { refName: 'GRCh38#0#chr1', start: 400, name: '400:0' }
    const runs = joinPieces([
      run(g, 'HG002#1#chr1', 100, [2]),
      walkRow(
        'HG002#1#chr1',
        3,
        1000,
        forward(7, 8),
        { abridged: { steps: 4, bp: 200, detail } },
        g.nodes,
      ),
      walkRow('HG002#1#chr1', 6, 2000, forward(3), {}, g.nodes),
    ])
    const { gapsOf } = walkCut(runs, g.nodes, 'GRCh38#0#chr1', [])
    expect(segments(gapsOf())).toEqual([
      [
        0,
        200,
        2000,
        [
          ['unread', 800, 200, 400],
          ['abridged', 200],
          ['unread', 800, 600, 200],
        ],
      ],
    ])
  })

  test('a 0 bp step a record attaches at the stretch’s end is kept as a 0 bp segment', () => {
    const g = graph()
    g.placeNode(9, {
      rank: 1,
      refName: 'HG002#1#chr1',
      start: 3000,
      end: 3000,
      passes: 0,
    })
    const { gapsOf } = walkCut(
      joinPieces([
        run(g, 'HG002#1#chr1', 100, [2, 5]),
        run(g, 'HG002#1#chr1', 500, [9, 3], { fragStart: 500 }),
      ]),
      g.nodes,
      'GRCh38#0#chr1',
      [record({})],
    )
    expect(segments(gapsOf())).toEqual([
      [
        0,
        200,
        500,
        [
          ['steps', 100, [5]],
          ['own', 200, 200],
          ['steps', 0, [9]],
        ],
      ],
    ])
  })

  test('records of 0 bp mark where the walk was inside unread pieces, and make no stretch at a row’s edge', () => {
    const g = graph()
    const runs = joinPieces([
      run(g, 'HG002#1#chr1', 100, [2]),
      walkRow('HG002#1#chr1', 5, 1000, forward(3), {}, g.nodes),
    ])
    const { gapsOf } = walkCut(runs, g.nodes, 'GRCh38#0#chr1', [
      record({ ownBp: 0, missingBp: 0 }),
      record({
        hapStart: 600,
        ownBp: 0,
        missingBp: 0,
        leavePos: 300,
        rejoinPos: 400,
      }),
    ])
    expect(segments(gapsOf())).toEqual([
      [
        0,
        200,
        1000,
        [
          ['unread', 400, 200, 300],
          ['own', 0, 0],
          ['unread', 400, 400, 200],
        ],
      ],
    ])
  })

  test('a row’s start reaches back over a record rejoining at it', () => {
    const g = graph()
    const { gapsOf } = walkCut(
      joinPieces([run(g, 'HG002#1#chr1', 500, [3], { fragStart: 500 })]),
      g.nodes,
      'GRCh38#0#chr1',
      [record({})],
    )
    expect(segments(gapsOf())).toEqual([[1, 200, 500, [['own', 300, 200]]]])
  })

  test('several records in one stretch each keep their own sequence, with what lies between them unread', () => {
    const g = graph()
    const runs = joinPieces([
      run(g, 'HG002#1#chr1', 100, [2]),
      walkRow('HG002#1#chr1', 5, 2000, forward(3), {}, g.nodes),
    ])
    const { gapsOf } = walkCut(runs, g.nodes, 'GRCh38#0#chr1', [
      record({ ownBp: 100, missingBp: 0, rejoinPos: 300 }),
      record({
        hapStart: 1000,
        ownBp: 100,
        missingBp: 0,
        leavePos: 400,
        rejoinPos: 500,
      }),
    ])
    expect(segments(gapsOf())).toEqual([
      [
        0,
        200,
        2000,
        [
          ['own', 100, 0],
          ['unread', 700, 300, 400],
          ['own', 100, 0],
          ['unread', 900, 500, 200],
        ],
      ],
    ])
  })

  test('between two fragments of one run in a contained cut lie the steps the cut read', () => {
    const g = graph()
    const { gapsOf } = walkCut(
      joinPieces([run(g, 'HG002#1#chr1', 0, [2, 3, 1, 2, 4])]),
      g.nodes,
      'GRCh38#0#chr1',
      [],
      { contained: true },
    )
    expect(segments(gapsOf()).filter(([edge]) => edge === 0)).toEqual([
      [0, 200, 300, [['read', 100]]],
    ])
  })
})
