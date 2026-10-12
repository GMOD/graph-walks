import { TabixIndexedFile } from '@gmod/tabix'

import { graphTablesGFA } from './testing.ts'
import { RegionTooLargeError, panSNHaplotype } from './util.ts'
import { WalkReader } from './walkReader.ts'
import { walkNameFilter } from './walkRows.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

function file(path: string) {
  return new TabixIndexedFile({ path, tbiPath: `${path}.tbi` })
}

// `<prefix>.walks.bed.gz` and, for format 1, the node and link files beside it
function reader(prefix: string, chunkWithoutHeader?: number) {
  const walks = file(`${DATA}${prefix}.walks.bed.gz`)
  return new WalkReader(
    prefix.startsWith('walks2_')
      ? { walks }
      : {
          walks,
          nodes: file(`${DATA}${prefix}.nodes.bed.gz`),
          links: file(`${DATA}${prefix}.links.bed.gz`),
        },
    chunkWithoutHeader,
  )
}

async function gfa(
  walks: WalkReader,
  start: number,
  end: number,
  haplotypes?: string[],
  refName = 'GRCh38#0#chr1',
) {
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

// Four walks over a 3.6 kb reference in 1 kb chunks, built by gfa-to-tabix
// --walks --refs GRCh38 --chunk 1000 (0.4.0):
//   W GRCh38  0 chr1 >1>2>3>4>5>6
//   W HG002   1 chr1 >1>2>7>4>5>6   7 replaces 3
//   W HG002   2 chr1 >1>2>3>4>8>6   8 replaces 5
//   W HG00097 1 chr1 >1>2>3>4>5>6
test('the header names the reference, the haplotypes and the chunk size', async () => {
  const header = await reader('walks_header').header()
  expect(header.references).toEqual(['GRCh38'])
  expect(header.haplotypes).toEqual(['HG00097#1', 'HG002#1', 'HG002#2'])
  expect(header.chunk).toBe(1000)
})

test('a cut keeps the reference and the haplotypes asked for', async () => {
  const names = walkNames(
    await gfa(reader('walks_header'), 1000, 2600, ['HG002#2']),
  )
  expect(names).toEqual(['GRCh38#0#chr1', 'HG002#2#chr1'])
})

// An unchopped graph whose node 1 spans 0-2600, filed under chunk 0 alone:
//   W GRCh38 0 chr1 >1>2>4
//   W HG002  1 chr1 >1>3>4   3 replaces 2
test("a node longer than a chunk is read back as far as the header's maxnode reaches", async () => {
  const text = await gfa(reader('walks_maxnode'), 2100, 2500)
  expect(text).toMatch(/^S\t1\t/m)
  expect(text).toMatch(/^W\tGRCh38\t0\tchr1\t0\t3200\t>1>2$/m)
  expect(walkNames(text)).toEqual(['GRCh38#0#chr1', 'HG002#1#chr1'])
})

test('a header without maxnode reads from the chunk before the window', async () => {
  const getLines = vi.spyOn(TabixIndexedFile.prototype, 'getLines')
  await gfa(reader('walks_header'), 2100, 2600)
  expect(new Set(getLines.mock.calls.map(call => call[1]))).toEqual(
    new Set([1000]),
  )
  getLines.mockRestore()
})

// Twelve 500 bp reference nodes in 1 kb chunks, whose rows name the chunks
// their path goes on in: HG002#1 leaves chunk 1000 for 20 and 21, 2 kb each,
// and rejoins at 9, in chunk 4000; HG003#1's contig ends on 22
test('a walk the read cuts off is read on to the chunk its row names', async () => {
  const getLines = vi.spyOn(TabixIndexedFile.prototype, 'getLines')
  const { fragments } = await reader('walks_follow').cut('GRCh38#0#chr1', {
    start: 1000,
    end: 1500,
  })
  const reads = getLines.mock.calls.map(call => [call[1], call[2]])
  getLines.mockRestore()
  expect(new Set(reads.map(String))).toEqual(new Set(['0,1500', '2000,4001']))
  const hg002 = fragments.find(f => f.name === 'HG002#1#chr1')!
  expect([...hg002.ids]).toContain(9)
})

test('a cut on a signal already aborted reads no row', async () => {
  const walks = reader('walks_header')
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
    reader('walks_header').cut('GRCh38#0#chr1', { start: NaN, end: 2600 }),
  ).rejects.toThrow("is not a finite range to query the graph's index for")
})

test('a cut over its step budget says how far to zoom in', async () => {
  const error = await reader('walks_header')
    .cut('GRCh38#0#chr1', { start: 0, end: 3600 }, { stepBudget: 2 })
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(RegionTooLargeError)
  expect((error as RegionTooLargeError).regionTooLarge).toBe(true)
})

test('a format 1 walk file without its node file says so', async () => {
  const walks = new WalkReader({
    walks: file(`${DATA}walks1_header.GRCh38.walks.bed.gz`),
  })
  await expect(
    walks.cut('GRCh38#0#chr1', { start: 0, end: 600 }),
  ).rejects.toThrow(/needs its node file/)
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
    expect(graph.sequences.get(id)).toHaveLength(node.end - node.start)
  }
  expect(await walks.referenceSequences()).toEqual(['GRCh38#0#chrM'])
})

// Each walks_* fixture built again by gfa-to-tabix 0.11.0 as format 1
// (walks1_*) and as format 2 (walks2_*). A format 2 cut is the format 1 cut.
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
])('walks_$name in format 2', ({ name, windows }) => {
  const one = reader(`walks1_${name}.GRCh38`)
  const two = reader(`walks2_${name}.GRCh38`)

  test('names the same haplotypes', async () => {
    expect((await two.header()).haplotypes).toEqual(
      (await one.header()).haplotypes,
    )
    expect((await two.header()).format).toBe(2)
  })

  test.each(windows.map(([start, end]) => ({ start: start!, end: end! })))(
    'cuts $start-$end as format 1 does',
    async ({ start, end }) => {
      for (const haplotypes of [
        undefined,
        ['HG002#1'],
        ['HG002'],
        ['HG00097', 'HG003'],
      ]) {
        expect(await gfa(two, start, end, haplotypes)).toEqual(
          await gfa(one, start, end, haplotypes),
        )
      }
    },
  )
})
