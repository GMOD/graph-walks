import { readFileSync } from 'node:fs'

import { TabixIndexedFile } from '@gmod/tabix'

import { graphTablesGFA } from './testing.ts'
import { RegionTooLargeError, panSNHaplotype } from './util.ts'
import { WalkFormatError, WalkReader } from './walkReader.ts'
import { walkNameFilter } from './walkRows.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

function reader(name: string) {
  const path = `${DATA}${name}.walks.bed.gz`
  return new WalkReader(new TabixIndexedFile({ path, tbiPath: `${path}.tbi` }))
}

async function gfa(
  walks: WalkReader,
  start: number,
  end: number,
  haplotypes?: string[],
) {
  const refName = 'GRCh38#0#chr1'
  const { graph, kept, fragments } = await walks.cut(
    refName,
    { start, end },
    { keep: walkNameFilter(haplotypes, panSNHaplotype(refName)) },
  )
  return graphTablesGFA(graph.tables(kept, fragments))
}

const walkNames = (text: string) =>
  text
    .split('\n')
    .filter(l => l.startsWith('W\t'))
    .map(l => l.split('\t').slice(1, 4).join('#'))

// Four walks over a 3.6 kb reference in 1 kb chunks:
//   W GRCh38  0 chr1 >1>2>3>4>5>6
//   W HG002   1 chr1 >1>2>7>4>5>6   7 replaces 3
//   W HG002   2 chr1 >1>2>3>4>8>6   8 replaces 5
//   W HG00097 1 chr1 >1>2>3>4>5>6
test('the header names the format, the reference, the haplotypes and the chunk size', async () => {
  const header = await reader('header.GRCh38').header()
  expect(header).toMatchObject({
    format: 5,
    chunk: 1000,
    references: ['GRCh38'],
    haplotypes: ['HG00097#1', 'HG002#1', 'HG002#2'],
  })
})

test('a format 1 walk file is refused, with what to rebuild it with', async () => {
  const walks = reader('format1')
  await expect(walks.header()).rejects.toThrow(WalkFormatError)
  await expect(
    walks.cut('GRCh38#0#chr1', { start: 0, end: 600 }),
  ).rejects.toThrow(/rebuild it with gfa-to-tabix 1.4/)
})

test.each([2, 3, 4])(
  'a format %i walk file is refused, with what to rebuild it with',
  async format => {
    const walks = new WalkReader({
      getHeaderLines: () =>
        Promise.resolve([`#walks\tformat:i:${format}\tchunk:i:1000`]),
      getLines: () => Promise.resolve(),
      getReferenceSequenceNames: () => Promise.resolve([]),
    })
    await expect(walks.header()).rejects.toThrow(
      `This walk file is format ${format}; @jbrowse/graph-walks reads format 5, so rebuild it with gfa-to-tabix 1.4 or later`,
    )
  },
)

test('a cut keeps the reference and the haplotypes asked for', async () => {
  const names = walkNames(
    await gfa(reader('header.GRCh38'), 1000, 2600, ['HG002#2']),
  )
  expect(names).toEqual(['GRCh38#0#chr1', 'HG002#2#chr1'])
})

// An unchopped graph whose node 1 spans 0-2600, filed under chunk 0 alone:
//   W GRCh38 0 chr1 >1>2>4
//   W HG002  1 chr1 >1>3>4   3 replaces 2
test("a node longer than a chunk is read back as far as the header's maxnode reaches", async () => {
  const text = await gfa(reader('maxnode.GRCh38'), 2100, 2500)
  expect(text).toMatch(/^S\t1\t/m)
  expect(text).toMatch(/^W\tGRCh38\t0\tchr1\t0\t3800\t>1>2>4$/m)
  expect(walkNames(text)).toEqual(['GRCh38#0#chr1', 'HG002#1#chr1'])
})

test('a cut reads from the chunk before the window', async () => {
  const getLines = vi.spyOn(TabixIndexedFile.prototype, 'getLines')
  await reader('header.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 2100, end: 2600 },
    { context: 0 },
  )
  expect(getLines.mock.calls[0]?.slice(1, 3)).toEqual([1000, 2600])
  getLines.mockRestore()
})

// Twelve 500 bp reference nodes in 1 kb chunks: HG002#1 leaves chunk 1000 for
// 20 and 21, 2 kb each, and rejoins at 9, in chunk 4000; HG003#1's contig
// ends on 22
test('a walk the read cuts off is read on to the chunk its row names', async () => {
  const getLines = vi.spyOn(TabixIndexedFile.prototype, 'getLines')
  const { fragments } = await reader('follow.GRCh38').cut('GRCh38#0#chr1', {
    start: 1000,
    end: 1500,
  })
  const reads = getLines.mock.calls.map(call => [call[1], call[2]])
  getLines.mockRestore()
  expect(new Set(reads.map(String))).toEqual(new Set(['0,1500', '2000,4001']))
  const hg002 = fragments.find(f => f.name === 'HG002#1#chr1')!
  expect([...hg002.ids]).toContain(9)
})

// Twenty 500 bp reference nodes in 1 kb chunks, built by gfa-to-tabix:
//   W GRCh38 0 chr1 0    10000 >1>2>...>20
//   W HG002  1 chr1 0    2000  >1>2>3>4
//   W HG002  1 chr1 2500 3200  >30          visits no reference node
//   W HG002  1 chr1 4000 7000  >7>8>...>12
//   W HG003  1 chr1 0    10000 >1>2>3>19>20>4>5>...>18
const records = async (start: number, end: number, name: string) => {
  const { fragments } = await reader('records.GRCh38').cut(
    'GRCh38#0#chr1',
    { start, end },
    { context: 0 },
  )
  return fragments
    .filter(f => f.name === name)
    .map(f => [f.fragStart, f.hapStart, f.hapEnd, f.before, f.after])
}
const chunk = (start: number) => ({ refName: 'GRCh38#0#chr1', start })

test("a fragment stops at its record's end with the bp to the contig's next record", async () => {
  expect(await records(1000, 2000, 'HG002#1#chr1')).toEqual([
    [0, 1000, 2000, { kind: 'end', trimmed: 1000 }, { kind: 'break', bp: 500 }],
  ])
  expect(await records(3000, 3600, 'HG002#1#chr1')).toEqual([
    [
      4000,
      4000,
      5000,
      { kind: 'break', bp: 800 },
      { kind: 'elsewhere', chunk: chunk(5000), trimmed: 1000 },
    ],
  ])
  expect(await records(5000, 6000, 'HG002#1#chr1')).toEqual([
    [
      4000,
      6000,
      7000,
      { kind: 'elsewhere', chunk: chunk(3000), trimmed: 1000 },
      { kind: 'end' },
    ],
  ])
})

test('a walk leaving the window for eight chunks stops where it leaves and comes back', async () => {
  expect(await records(1000, 2000, 'HG003#1#chr1')).toEqual([
    [
      0,
      1000,
      1500,
      { kind: 'end', trimmed: 1000 },
      { kind: 'elsewhere', chunk: chunk(9000), bp: 1000 },
    ],
    [
      0,
      2500,
      3000,
      { kind: 'elsewhere', chunk: chunk(9000), bp: 1000 },
      { kind: 'elsewhere', chunk: chunk(3000), trimmed: 1000 },
    ],
  ])
})

test('a cut on a signal already aborted reads no row', async () => {
  const walks = reader('header.GRCh38')
  await walks.header()
  const getLines = vi.spyOn(TabixIndexedFile.prototype, 'getLines')
  const aborted = new AbortController()
  aborted.abort(new Error('panned away'))
  await expect(
    walks.cut(
      'GRCh38#0#chr1',
      { start: 1000, end: 2600 },
      { signal: aborted.signal },
    ),
  ).rejects.toThrow('panned away')
  expect(getLines).not.toHaveBeenCalled()
  getLines.mockRestore()
})

test('a range that is not finite queries nothing', async () => {
  await expect(
    reader('header.GRCh38').cut('GRCh38#0#chr1', { start: NaN, end: 2600 }),
  ).rejects.toThrow("is not a finite range to query the graph's index for")
})

test('a cut over its step budget says how far to zoom in', async () => {
  const error = await reader('header.GRCh38')
    .cut('GRCh38#0#chr1', { start: 0, end: 3600 }, { stepBudget: 2 })
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(RegionTooLargeError)
  expect((error as RegionTooLargeError).regionTooLarge).toBe(true)
})

// The chrM graph of HPRC's release 2, from `vg convert -f` of its GBZ and
// gfa-to-tabix --walks --refs GRCh38 --sequences
test("a cut keeps each node's bases when asked", async () => {
  const walks = reader('hprc-chrM.GRCh38')
  const { graph, kept, fragments } = await walks.cut(
    'GRCh38#0#chrM',
    { start: 8000, end: 8300 },
    { sequences: true },
  )
  expect(walkNames(graphTablesGFA(graph.tables(kept, fragments)))).toHaveLength(
    44,
  )
  for (const id of kept) {
    const node = graph.nodes.get(id)!
    expect(graph.sequences.get(id)).toMatch(/^[ACGTN]+$/)
    expect(graph.sequences.get(id)).toHaveLength(node.end - node.start)
  }
  expect(await walks.referenceSequences()).toEqual(['GRCh38#0#chrM'])
})

test('a contig names its sequence on the reference sample', async () => {
  const walks = reader('hprc-chrM.GRCh38')
  expect(await walks.sequenceFor('GRCh38', 'chrM')).toBe('GRCh38#0#chrM')
  expect(await walks.sequenceFor('hg38', 'GRCh38#0#chrM')).toBe('GRCh38#0#chrM')
  expect(await walks.sequenceFor('GRCh38', 'chr1')).toBeUndefined()
})

test("a region cut names the file's sequence for a contig, and keeps the haplotypes asked for", async () => {
  const cut = await reader('header.GRCh38').cutRegion(
    { refName: 'chr1', start: 1000, end: 2600 },
    { haplotypes: ['HG002#2'] },
  )
  expect(cut.refName).toBe('GRCh38#0#chr1')
  expect(cut.sample).toBe('GRCh38')
  expect(walkNames(graphTablesGFA(cut.tables))).toEqual([
    'GRCh38#0#chr1',
    'HG002#2#chr1',
  ])
  const whole = await reader('header.GRCh38').cutRegion({
    refName: 'GRCh38#0#chr1',
    start: 1000,
    end: 2600,
  })
  expect(graphTablesGFA(whole.tables)).toBe(
    await gfa(reader('header.GRCh38'), 1000, 2600),
  )
})

test('a region cut on a sequence the file lacks names the ones it holds', async () => {
  await expect(
    reader('header.GRCh38').cutRegion({ refName: 'chr2', start: 0, end: 100 }),
  ).rejects.toThrow(
    'The walk file has no sequence chr2 on GRCh38; it holds GRCh38#0#chr1',
  )
})

test('a region cut over its node limit says how far to zoom in', async () => {
  const error = await reader('header.GRCh38')
    .cutRegion({ refName: 'chr1', start: 0, end: 3600 }, { nodeLimit: 4 })
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(RegionTooLargeError)
  expect(error).toMatchObject({ limit: 4, windowBp: 3600 })
  expect((error as Error).message).toMatch(
    /^Zoom in to about [\d.]+K?bp, or choose fewer haplotypes, to see the graph$/,
  )
})

// A walk file from rows written space-separated, served in chunk order or
// reversed
function inline(lines: string[], reversed = false) {
  const rows = lines.map(line => line.split(' ').join('\t'))
  const header = rows.filter(row => row.startsWith('#'))
  const body = rows.filter(row => !row.startsWith('#'))
  return new WalkReader({
    getHeaderLines: () => Promise.resolve(header),
    getLines: (refName, start, end, { lineCallback }) => {
      const hits = body.filter(row => {
        const [seq, cs] = row.split('\t')
        return seq === refName && +cs! < end && +cs! + 1 > start
      })
      for (const row of reversed ? hits.reverse() : hits) {
        lineCallback(row)
      }
      return Promise.resolve()
    },
    getReferenceSequenceNames: () => Promise.resolve(['R#0#chr']),
  })
}

// Eight 4 bp reference nodes in 16 bp chunks:
//   W R  0 chr >1>2>3>4>5>6>7>8
//   W H1 1 c   >1>2>6>3>4>5>6>7>8
//   W H2 1 c   >5>6>2>7>8
// `pieces` files each visit to the other chunk as a piece of its own. `stays`
// keeps the visit in its piece as a route interior, as gfa-to-tabix's --stay
// does on a reference that passes 2 and 6 twice, so chunk 0's route steps on
// 6, which chunk 16's ref row places, and chunk 16's on 2.
const tinyHeader = [
  '#walks format:i:5 chunk:i:16 maxnode:i:4',
  '#reference R',
  '#haplotype H1#1',
  '#haplotype H2#1',
]
const pieces = [
  ...tinyHeader,
  'R#0#chr 0 1 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 nx:i:16',
  'R#0#chr 0 1 walk H1#1#c 0 0 0 2 2 2 . nx:i:16',
  'R#0#chr 0 1 walk H1#1#c 0 12 2 2 6 2 . pv:i:16 nx:i:16',
  'R#0#chr 0 1 walk H2#1#c 0 8 1 1 4 0 . pv:i:16 nx:i:16',
  'R#0#chr 16 17 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 pv:i:0',
  'R#0#chr 16 17 walk H1#1#c 0 8 1 1 12 0 . pv:i:0 nx:i:0',
  'R#0#chr 16 17 walk H1#1#c 0 20 3 4 10 6 . pv:i:0',
  'R#0#chr 16 17 walk H2#1#c 0 0 0 2 10 2 . nx:i:0',
  'R#0#chr 16 17 walk H2#1#c 0 12 2 2 14 2 . pv:i:0',
]
const stays = [
  ...tinyHeader,
  'R#0#chr 0 1 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 nx:i:16',
  'R#0#chr 0 1 nodes 6 4',
  'R#0#chr 0 1 route 4 8 -6',
  'R#0#chr 0 1 walk H1#1#c 0 0 0 5 2 6 0 nx:i:16',
  'R#0#chr 16 17 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 pv:i:0',
  'R#0#chr 16 17 nodes 2 4',
  'R#0#chr 16 17 route 12 -8 10',
  'R#0#chr 16 17 walk H1#1#c 0 20 1 4 10 6 . pv:i:0',
  'R#0#chr 16 17 walk H2#1#c 0 0 0 5 10 6 0',
]

// Each window reads both chunks, the first as the chunk before it
describe.each([
  [4, 20],
  [16, 32],
  [20, 28],
  [0, 32],
])('a reference node a route steps on, cut %i-%i', (start, end) => {
  test.each([
    { haplotypes: undefined },
    { haplotypes: ['H1'] },
    { haplotypes: ['H2'] },
  ])(
    'keeps its reference position, haplotypes $haplotypes',
    async ({ haplotypes }) => {
      const cut = (lines: string[], reversed?: boolean) =>
        inline(lines, reversed)
          .cutRegion({ refName: 'chr', start, end }, { haplotypes, context: 0 })
          .then(({ tables }) => graphTablesGFA(tables))
      const text = await cut(pieces)
      expect(text).toMatch(/^S\t2\t\*\tLN:i:4\tSN:Z:R#0#chr\tSO:i:4\tSR:i:0$/m)
      expect(text).toMatch(/^S\t6\t\*\tLN:i:4\tSN:Z:R#0#chr\tSO:i:20\tSR:i:0$/m)
      expect(await cut(stays)).toBe(text)
      expect(await cut(stays, true)).toBe(text)
    },
  )
})

// The GFA format 1 files of the same graphs cut to, recorded before format 1
// was dropped, when a format 2 cut was checked to be the format 1 cut
const goldens = JSON.parse(
  readFileSync(`${DATA}goldens.json`, 'utf8'),
) as Record<string, string>

describe.each([
  {
    name: 'header',
    windows: [
      [0, 600],
      [1000, 2600],
      [2100, 3600],
    ],
  },
  {
    name: 'maxnode',
    windows: [
      [2100, 2500],
      [0, 3800],
    ],
  },
  {
    name: 'rejoin',
    windows: [
      [1000, 1500],
      [0, 6000],
      [4000, 4200],
    ],
  },
  {
    name: 'follow',
    windows: [
      [1000, 1500],
      [0, 6000],
      [3500, 5200],
    ],
  },
])('$name', ({ name, windows }) => {
  test.each(windows.map(([start, end]) => ({ start: start!, end: end! })))(
    'cuts $start-$end as recorded',
    async ({ start, end }) => {
      const walks = reader(`${name}.GRCh38`)
      for (const haplotypes of [
        undefined,
        ['HG002#1'],
        ['HG002'],
        ['HG00097', 'HG003'],
      ]) {
        const key = `${name} ${start}-${end} ${haplotypes?.join(',') ?? 'all'}`
        expect(await gfa(walks, start, end, haplotypes)).toBe(goldens[key])
      }
    },
  )
})

// A walk file from rows written space-separated, which notes the chunk of
// each row it serves
function rowsReader(lines: string[], served: number[] = []) {
  const rows = lines.map(line => line.split(' ').join('\t'))
  return new WalkReader({
    getHeaderLines: () =>
      Promise.resolve(rows.filter(row => row.startsWith('#'))),
    getLines: (refName, start, end, { lineCallback }) => {
      for (const row of rows) {
        const [seq, cs] = row.split('\t')
        if (seq === refName && +cs! < end && +cs! + 1 > start) {
          served.push(+cs!)
          lineCallback(row)
        }
      }
      return Promise.resolve()
    },
    getReferenceSequenceNames: () => Promise.resolve(['R#0#chr']),
  })
}

// Twelve 4 bp reference nodes in 16 bp chunks; H#1 goes from 7 to 11, in
// chunk 32, and back to 6:
//   W R  0 chr >1>2>3>4>5>6>7>8>9>10>11>12
//   W H  1 h   >5>6>7>11>6>7>8
const copyJump = [
  '#walks format:i:5 chunk:i:16 maxnode:i:4',
  '#reference R',
  '#haplotype H#1',
  'R#0#chr 0 1 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 nx:i:16',
  'R#0#chr 16 17 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 pv:i:0 nx:i:32',
  'R#0#chr 16 17 walk H#1#h 0 0 0 3 10 4 . nx:i:32',
  'R#0#chr 16 17 walk H#1#h 0 16 2 3 12 4 . pv:i:32',
  'R#0#chr 32 33 ref R#0#chr 0 32 2 4 18,2,2,2 4,4,4,4 pv:i:16',
  'R#0#chr 32 33 walk H#1#h 0 12 1 1 22 0 . pv:i:16 nx:i:16',
]
const copyJumpCut = async (
  start: number,
  end: number,
  contained: boolean,
  served: number[] = [],
) => {
  const { graph, kept, fragments } = await rowsReader(copyJump, served).cut(
    'R#0#chr',
    { start, end },
    { context: 0, contained, nodeLimit: 5 },
  )
  return {
    text: graphTablesGFA(graph.tables(kept, fragments)),
    walks: fragments
      .filter(f => f.name === 'H#1#h')
      .map(f => [f.hapStart, f.hapEnd, [...f.ids], f.before, f.after]),
  }
}

test('a contained cut splits a walk where it steps onto the reference off the window', async () => {
  const detour = await copyJumpCut(20, 36, false)
  expect(detour.walks).toEqual([
    [4, 28, [6, 7, 11, 6, 7, 8], { kind: 'end', trimmed: 4 }, { kind: 'end' }],
  ])
  expect(detour.text).toMatch(/^S\t11\t/m)
  const contained = await copyJumpCut(20, 36, true)
  expect(contained.walks).toEqual([
    [4, 12, [6, 7], { kind: 'end', trimmed: 4 }, { kind: 'end', trimmed: 16 }],
    [16, 28, [6, 7, 8], { kind: 'end', trimmed: 16 }, { kind: 'end' }],
  ])
  expect(contained.text).not.toMatch(/^S\t11\t/m)
  expect(contained.text).toMatch(/^W\tR\t0\tchr\t20\t36\t>6>7>8>9$/m)
})

test('a contained cut reads on only to chunks under the window', async () => {
  const served: number[] = []
  const { walks } = await copyJumpCut(20, 28, true, served)
  const elsewhere = {
    kind: 'elsewhere',
    chunk: { refName: 'R#0#chr', start: 32 },
    bp: 4,
  }
  expect(walks).toEqual([
    [4, 12, [6, 7], { kind: 'end', trimmed: 4 }, elsewhere],
    [16, 24, [6, 7], elsewhere, { kind: 'end', trimmed: 4 }],
  ])
  expect(new Set(served)).toEqual(new Set([0, 16]))
  const detour: number[] = []
  await copyJumpCut(20, 28, false, detour)
  expect(new Set(detour)).toEqual(new Set([0, 16, 32]))
})

// Eight chunks of four 4 bp reference nodes
const eightChunks = [
  '#walks format:i:5 chunk:i:16 maxnode:i:4',
  '#reference R',
  ...Array.from({ length: 8 }, (_, c) =>
    [
      `R#0#chr ${16 * c} ${16 * c + 1} ref R#0#chr 0 ${16 * c} ${c} 4`,
      `${8 * c + 2},2,2,2 4,4,4,4`,
      ...(c > 0 ? [`pv:i:${16 * (c - 1)}`] : []),
      ...(c < 7 ? [`nx:i:${16 * (c + 1)}`] : []),
    ].join(' '),
  ),
]

test('a cut over its node limit stops reading at the chunk that puts it over', async () => {
  const served: number[] = []
  const error = await rowsReader(eightChunks, served)
    .cutRegion(
      { refName: 'chr', start: 0, end: 128 },
      { context: 0, nodeLimit: 6 },
    )
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(RegionTooLargeError)
  expect(error).toMatchObject({ limit: 6, windowBp: 128, fitsBp: 24 })
  expect(served).toEqual([0, 16, 32])
  const whole = await rowsReader(eightChunks).cutRegion(
    { refName: 'chr', start: 0, end: 128 },
    { context: 0, nodeLimit: 32 },
  )
  expect(whole.tables.nodes.names).toHaveLength(32)
})

// gfa-to-tabix's walks-two-refs.gfa on R in 10 bp chunks, once as written and
// once with --abridge-bytes 1, which abridges every walk in every chunk
async function twoRefs(name: string, expand?: (name: string) => boolean) {
  const { graph, kept, fragments } = await reader(name).cut(
    'R#0#chr',
    { start: 0, end: 30 },
    { expand },
  )
  return { tables: graph.tables(kept, fragments), fragments }
}

test('an abridged walk cuts to a fragment from where it leaves R to where it rejoins, with its steps and bp', async () => {
  const { tables, fragments } = await twoRefs('abridged.R')
  const h = fragments.filter(f => f.name === 'H#1#h')
  expect(h.map(f => [f.hapStart, f.hapEnd, f.abridged?.steps])).toEqual([
    [100, 112, 2],
    [112, 127, 3],
    [127, 133, 1],
  ])
  const { names, abridgedSteps } = tables.walks
  names.forEach((name, i) => {
    expect(Number.isNaN(abridgedSteps![i])).toBe(name === 'R#0#chr')
  })
})

test('an expanded walk is read from its detail records as the unabridged file cuts it', async () => {
  const all = () => true
  expect(graphTablesGFA((await twoRefs('abridged.R', all)).tables)).toBe(
    graphTablesGFA((await twoRefs('unabridged.R')).tables),
  )
})

test('the reference sequences leave out detail records', async () => {
  expect(await reader('abridged.R').referenceSequences()).toEqual(['R#0#chr'])
})
