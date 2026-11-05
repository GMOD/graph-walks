import { graphTablesGFA } from './testing.ts'
import {
  WalkGraph,
  byteBudgetError,
  chunkQueryStart,
  joinPieces,
  lookbackChunks,
  stepBudgetError,
  walkCut,
  walkHeader,
  walkNameFilter,
  walkQueryStart,
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

test('a cut queries from the start of the chunk before the window', () => {
  expect(chunkQueryStart(20_000_000, 65_536)).toBe(19_922_944)
  expect(chunkQueryStart(65_536, 65_536)).toBe(0)
  expect(chunkQueryStart(10, 65_536)).toBe(0)
})

test("a walk file's header states its format, chunk, references, haplotypes and holes", () => {
  expect(walkHeader(['#walks\tchunk:i:65536'])).toEqual({
    chunk: 65_536,
    references: [],
    haplotypes: [],
    long: [],
    holes: [],
  })
  expect(walkHeader(['#walks']).chunk).toBeUndefined()
  expect(
    walkHeader([
      '#walks\tformat:i:6\tchunk:i:1000\tmaxnode:i:1024\tstay:i:311',
      '#reference\tGRCh38',
      '#reference\tCHM13',
      '#haplotype\tHG00097#1',
      '#haplotype\tHG002#1',
      '#hole\tGRCh38#0#chr1\t12\t20',
    ]),
  ).toEqual({
    format: 6,
    chunk: 1000,
    maxNode: 1024,
    references: ['GRCh38', 'CHM13'],
    haplotypes: ['HG00097#1', 'HG002#1'],
    long: [],
    holes: [{ refName: 'GRCh38#0#chr1', start: 12, end: 20 }],
  })
})

test('a cut looks back as many chunks as the longest node spans', () => {
  expect(lookbackChunks(undefined, 65_536)).toBe(1)
  expect(lookbackChunks(0, 65_536)).toBe(1)
  expect(lookbackChunks(1024, 65_536)).toBe(1)
  expect(lookbackChunks(65_536, 65_536)).toBe(1)
  expect(lookbackChunks(65_537, 65_536)).toBe(2)
  expect(lookbackChunks(300_000, 65_536)).toBe(5)
  expect(chunkQueryStart(20_000_000, 65_536, 5)).toBe(19_660_800)
  expect(chunkQueryStart(200_000, 65_536, 5)).toBe(0)
})

// GRCh38#0#chr22 has a 1 Mb N run at 10,510,000; a window reads past the
// chunk before its own only where it overlaps that run
test('a cut reads back to a long node only where the window overlaps it', () => {
  const header = {
    ...walkHeader([
      '#walks\tformat:i:6\tchunk:i:65536\tmaxnode:i:1000000',
      '#long\tGRCh38#0#chr22\t10510000\t11510000',
    ]),
    chunk: 65_536,
  }
  expect(header.long).toEqual([
    { refName: 'GRCh38#0#chr22', start: 10_510_000, end: 11_510_000 },
  ])
  const at = (refName: string, start: number) =>
    walkQueryStart(header, refName, start, start + 1000)
  expect(at('GRCh38#0#chr22', 20_000_000)).toBe(19_922_944)
  expect(at('GRCh38#0#chr22', 11_000_000)).toBe(10_420_224)
  expect(at('GRCh38#0#chr21', 11_000_000)).toBe(10_878_976)
  expect(
    walkQueryStart(
      { ...walkHeader(['#walks\tmaxnode:i:1000000']), chunk: 65_536 },
      'GRCh38#0#chr22',
      20_000_000,
      20_001_000,
    ),
  ).toBe(18_939_904)
})

describe('the byte budget', () => {
  const chunk = 100
  // bytes a read from the query start to `end` fetches, 10 per chunk row
  // filed under the chunk's first base, and 80 under the chunk at 300
  const bytesTo = (end: number) => {
    let bytes = 0
    for (let cs = 0; cs < end; cs += chunk) {
      bytes += cs === 300 ? 80 : 10
    }
    return Promise.resolve(bytes)
  }

  test('passes a window under it', async () => {
    expect(
      await byteBudgetError(bytesTo, 40, { start: 150, end: 250 }, chunk),
    ).toBeUndefined()
  })

  test('names the span from the window start that fits', async () => {
    const error = (await byteBudgetError(
      bytesTo,
      50,
      { start: 150, end: 450 },
      chunk,
    ))!
    expect(error.name).toBe('NodeLimitError')
    expect(error.regionTooLarge).toBe(true)
    expect(error.fitsBp).toBe(150)
    expect(error.zoomHelps).toBe(true)
    expect(error.bytes).toBe(120)
    expect(error.message).toBe('Zoom in to about 150bp to see the graph')
  })

  test('says so where no zoom fits', async () => {
    const error = (await byteBudgetError(
      async end => (await bytesTo(end)) * 100_000,
      5_000_000,
      { start: 310, end: 320 },
      chunk,
    ))!
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

test("a record's own steps come from its leave piece, and a leave piece the cut lacks is asked for", () => {
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
  expect(lacking.unattached).toEqual([{ departure: record(200), base: 199 }])
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
