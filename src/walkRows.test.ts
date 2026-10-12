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
} from './walkRows.ts'

import type { WalkRow } from './walkRows.ts'

// A decoded walk row: a piece of `name` from `hapOffset`, through `steps`
function walkRow(
  name: string,
  piece: number,
  hapOffset: number,
  steps: [number, boolean][],
  more: Partial<WalkRow> = {},
): WalkRow {
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
    },
    ...more,
  }
}

const forward = (...ids: number[]) =>
  ids.map(id => [id, false] as [number, boolean])

test('a run takes the chunks its first and last rows say its path goes on in', () => {
  const rows = [
    walkRow('HG002#1#chr1', 0, 0, forward(10, 11), {
      next: { refName: 'GRCh38#0#chr1', start: 1000 },
    }),
    walkRow('HG002#1#chr1', 1, 300, forward(12), {
      prev: { refName: 'GRCh38#0#chr1', start: 0 },
      next: { refName: 'GRCh38#0#chr2', start: 5000 },
    }),
  ]
  const [run] = joinPieces(rows)
  expect([run!.prev, run!.next]).toEqual([
    undefined,
    { refName: 'GRCh38#0#chr2', start: 5000 },
  ])
})

test('continuation rows join by piece index, and a missing piece splits the walk', () => {
  const rows = [
    walkRow('HG002#1#chr1', 3, 900, forward(30)),
    walkRow('HG002#1#chr1', 1, 300, forward(12, 13)),
    walkRow('HG002#1#chr1', 0, 100, forward(10, 11)),
  ]
  const runs = joinPieces(rows)
  expect(runs.map(r => [r.hapOffset, [...r.ids]])).toEqual([
    [100, [10, 11, 12, 13]],
    [900, [30]],
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

describe('the step budget', () => {
  const chunk = 100
  const rows = (perChunk: [number, number][]): WalkRow[] =>
    perChunk.map(([chunkStart, n]) => ({
      name: 'HG002#1#chr1',
      fragStart: 0,
      hapOffset: 0,
      piece: 0,
      n,
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
      'Too many haplotype steps here to draw these haplotypes (50 steps read, over 2 × the step budget of 5); choose fewer haplotypes',
    )
  })
})

test('a cut queries from the start of the chunk before the window', () => {
  expect(chunkQueryStart(20_000_000, 65_536)).toBe(19_922_944)
  expect(chunkQueryStart(65_536, 65_536)).toBe(0)
  expect(chunkQueryStart(10, 65_536)).toBe(0)
})

test("a walk file's header states its chunk, references and haplotypes", () => {
  expect(walkHeader(['#walks\tchunk:i:65536'])).toEqual({
    chunk: 65_536,
    references: [],
    haplotypes: [],
  })
  expect(walkHeader(['#walks']).chunk).toBeUndefined()
  expect(
    walkHeader([
      '#walks\tchunk:i:1000',
      '#reference\tGRCh38',
      '#reference\tCHM13',
      '#haplotype\tHG00097#1',
      '#haplotype\tHG002#1',
    ]),
  ).toEqual({
    chunk: 1000,
    references: ['GRCh38', 'CHM13'],
    haplotypes: ['HG00097#1', 'HG002#1'],
  })
})

test("a 0.5.0 header's first line gives the longest node beside the chunk", () => {
  expect(
    walkHeader([
      '#walks\tchunk:i:65536\tmaxnode:i:1024\tcap:i:8192',
      '#reference\tGRCh38',
    ]),
  ).toEqual({
    chunk: 65_536,
    maxNode: 1024,
    references: ['GRCh38'],
    haplotypes: [],
  })
  expect(walkHeader(['#walks\tmaxnode:i:300000\tchunk:i:1000'])).toMatchObject({
    chunk: 1000,
    maxNode: 300_000,
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
    expect(error.message).toBe('Zoom in to about 150bp to see the graph')
  })

  test('says so where no zoom fits', async () => {
    const error = (await byteBudgetError(
      async end => (await bytesTo(end)) * 100_000,
      5_000_000,
      { start: 310, end: 320 },
      chunk,
    ))!
    expect(error.message).toBe(
      'Too much graph here to fetch (11.0 MB against a byte budget of 5.0 MB)',
    )
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
    g.placeNode(id, { rank, refName, start, end: start + 100 })
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

test('a node is placed once', () => {
  const g = graph()
  expect(g.nodes.size).toBe(6)
  expect(g.nodes.get(2)).toEqual({
    rank: 0,
    refName: 'GRCh38#0#chr1',
    start: 100,
    end: 200,
    onWindow: true,
  })
})

test('a cut is written as the GFA its tables stand for', () => {
  const g = graph()
  const runs = joinPieces([
    walkRow('GRCh38#0#chr1', 0, 0, forward(1, 2, 3, 4)),
    // brings allele 5 into the cut between two steps on the window
    walkRow('HG002#1#chr1', 0, 0, forward(2, 5, 3)),
    // enters at 3 from 1 through allele 5
    walkRow('HG002#2#chr1', 0, 1000, forward(6, 1, 5, 3, 4)),
  ])
  const { kept, fragments } = walkCut(runs, g.nodes, 'GRCh38#0#chr1')
  expect([...kept].sort()).toEqual([1, 2, 3, 5])
  expect(
    fragments.map(f => [f.name, f.hapStart, f.hapEnd, [...f.ids]]),
  ).toEqual([
    ['GRCh38#0#chr1', 0, 300, [1, 2, 3]],
    ['HG002#1#chr1', 0, 300, [2, 5, 3]],
    ['HG002#2#chr1', 1100, 1400, [1, 5, 3]],
  ])
  expect(graphTablesGFA(g.tables(kept, fragments)).split('\n')).toEqual([
    'H\tVN:Z:1.1',
    'S\t1\t*\tLN:i:100\tSN:Z:GRCh38#0#chr1\tSO:i:0\tSR:i:0',
    'S\t2\t*\tLN:i:100\tSN:Z:GRCh38#0#chr1\tSO:i:100\tSR:i:0',
    'S\t3\t*\tLN:i:100\tSN:Z:GRCh38#0#chr1\tSO:i:200\tSR:i:0',
    'S\t5\t*\tLN:i:100\tSN:Z:HG002#1#chr1\tSO:i:1000\tSR:i:1',
    'L\t1\t+\t2\t+\t0M',
    'L\t1\t+\t5\t+\t0M',
    'L\t2\t+\t3\t+\t0M',
    'L\t2\t+\t5\t+\t0M',
    'L\t5\t+\t3\t+\t0M',
    'W\tGRCh38\t0\tchr1\t0\t300\t>1>2>3',
    'W\tHG002\t1\tchr1\t0\t300\t>2>5>3',
    'W\tHG002\t2\tchr1\t1100\t1400\t>1>5>3',
  ])
})

test('a walk leaving the window is followed to where it rejoins the reference', () => {
  const g = graph()
  const runs = joinPieces([
    walkRow('GRCh38#0#chr1', 0, 0, forward(1, 2, 3, 4)),
    // leaves at 2 through allele 5 and rejoins at 4
    walkRow('HG002#1#chr1', 0, 0, forward(1, 2, 5, 4)),
    // deletes 3, from 2 straight to 4
    walkRow('HG002#2#chr1', 0, 0, forward(1, 2, 4)),
    // enters at 3 from 1 through allele 6, and runs on back into allele 5,
    // which HG002#1 brings into the cut
    walkRow('HG002#3#chr1', 0, 0, forward(5, 1, 6, 3)),
    // leaves at 2 through allele 5 and ends there
    walkRow('HG002#4#chr1', 0, 0, forward(2, 5)),
    // starts on allele 6 and enters at 3
    walkRow('HG002#5#chr1', 0, 0, forward(6, 3)),
  ])
  const { kept, fragments } = walkCut(runs, g.nodes, 'GRCh38#0#chr1')
  expect([...kept].sort()).toEqual([1, 2, 3, 4, 5, 6])
  expect(
    fragments.map(f => [f.name, f.hapStart, f.hapEnd, [...f.ids]]),
  ).toEqual([
    ['GRCh38#0#chr1', 0, 400, [1, 2, 3, 4]],
    ['HG002#1#chr1', 0, 400, [1, 2, 5, 4]],
    ['HG002#2#chr1', 0, 300, [1, 2, 4]],
    ['HG002#3#chr1', 0, 400, [5, 1, 6, 3]],
    ['HG002#4#chr1', 0, 200, [2, 5]],
    ['HG002#5#chr1', 0, 200, [6, 3]],
  ])
})

test('a walk is not followed to another reference sequence', () => {
  const g = graph()
  g.placeNode(7, { rank: 0, refName: 'GRCh38#0#chr2', start: 0, end: 100 })
  const runs = joinPieces([walkRow('HG002#1#chr1', 0, 0, forward(2, 3, 5, 7))])
  const { fragments } = walkCut(runs, g.nodes, 'GRCh38#0#chr1')
  expect(fragments.map(f => [...f.ids])).toEqual([[2, 3, 5]])
})

test('the reference walk is written first, whatever order the rows came in', () => {
  const g = graph()
  const runs = joinPieces([
    walkRow('HG002#2#chr1', 0, 0, forward(2, 3)),
    walkRow('GRCh38#0#chr1', 0, 0, forward(1, 2, 3, 4)),
    walkRow('HG002#1#chr1', 0, 0, forward(2, 5, 3)),
  ])
  const { kept, fragments } = walkCut(runs, g.nodes, 'GRCh38#0#chr1')
  expect(g.tables(kept, fragments).walks.names).toEqual([
    'GRCh38#0#chr1',
    'HG002#1#chr1',
    'HG002#2#chr1',
  ])
})
