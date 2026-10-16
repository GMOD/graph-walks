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
    format: 3,
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
  ).rejects.toThrow(/rebuild it with gfa-to-tabix 1.1/)
})

test.each([2, 4])(
  'a format %i walk file is refused, with what to rebuild it with',
  async format => {
    const walks = new WalkReader({
      getHeaderLines: () =>
        Promise.resolve([`#walks\tformat:i:${format}\tchunk:i:1000`]),
      getLines: () => Promise.resolve(),
      getReferenceSequenceNames: () => Promise.resolve([]),
    })
    await expect(walks.header()).rejects.toThrow(
      `This walk file is format ${format}; @jbrowse/graph-walks reads format 3, so rebuild it with gfa-to-tabix 1.1 or later`,
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
