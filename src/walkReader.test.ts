import { BlobFile, LocalFile } from 'generic-filehandle2'

import { bgzf, graphTablesGFA, walkFileOf } from './testing.ts'
import { RegionTooLargeError, panSNHaplotype } from './util.ts'
import { WalkRowError } from './walkChunk.ts'
import { departureKind, placementChoice } from './walkDeparture.ts'
import { WalkFormatError } from './walkFile.ts'
import {
  DepartureReadError,
  WALK_BUDGETS,
  WalkReader,
  readInOrder,
} from './walkReader.ts'
import { walkNameFilter } from './walkRows.ts'

import type { WalkDeparture } from './walkDeparture.ts'
import type { WalkFragment } from './walkRows.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

function reader(name: string) {
  const path = `${DATA}${name}.walks.gz`
  return WalkReader.fromFilehandles({
    data: new LocalFile(path),
    directory: new LocalFile(`${path}.idx`),
  })
}

// The units of `name`'s sequence each read in `run` fetches, as
// `kind name`, and how many reads it makes
async function unitsRead(
  name: string,
  refName: string,
  run: (walks: WalkReader) => Promise<unknown>,
) {
  const walks = reader(name)
  const units = await walks.units(refName)
  await walks.header()
  const read = vi.spyOn(LocalFile.prototype, 'read')
  try {
    await run(walks)
    return read.mock.calls.map(([length, position = 0]) =>
      units
        .filter(u => u.offset >= position && u.offset < position + length)
        .map(u => `${u.kind} ${u.name}`),
    )
  } finally {
    read.mockRestore()
  }
}

async function gfa(
  walks: WalkReader,
  start: number,
  end: number,
  haplotypes?: string[],
) {
  const refName = 'GRCh38#0#chr1'
  const { graph, kept, fragments, departures } = await walks.cut(
    refName,
    { start, end },
    { keep: walkNameFilter(haplotypes, panSNHaplotype(refName)) },
  )
  return graphTablesGFA(graph.tables(kept, fragments, departures))
}

const walkNames = (text: string) =>
  text
    .split('\n')
    .filter(l => l.startsWith('W\t'))
    .map(l => l.split('\t').slice(1, 4).join('#'))

const fragment = (f: WalkFragment) => [
  f.name,
  f.hapStart,
  f.hapEnd,
  [...f.ids].map((id, k) => (f.rev[k] ? '<' : '>') + id).join(''),
  f.trimmedBefore,
  f.trimmedAfter,
]

const departure = (d: WalkDeparture) => [
  d.name,
  d.hapStart,
  d.ownBp,
  d.missingBp,
  `${d.leavePos}${d.leaveReversed ? '-' : '+'}`,
  `${d.rejoinPos}${d.rejoinReversed ? '-' : '+'}`,
  departureKind(d),
  d.fromSteps ? 'steps' : 'record',
  d.steps
    ? [...d.steps.ids]
        .map((id, k) => (d.steps!.rev[k] ? '<' : '>') + id)
        .join('')
    : undefined,
]

// Four walks over a 3.6 kb reference in 1 kb chunks (header.GRCh38.gfa):
//   W GRCh38  0 chr1 >1>2>3>4>5>6
//   W HG002   1 chr1 >1>2>7>4>5>6   7 replaces 3
//   W HG002   2 chr1 >1>2>3>4>8>6   8 replaces 5
//   W HG00097 1 chr1 >1>2>3>4>5>6
test('the header names the format, the reference, the haplotypes and the tile size', async () => {
  const header = await reader('header.GRCh38').header()
  expect(header).toMatchObject({
    format: 7,
    tile: 1000,
    references: ['GRCh38'],
    haplotypes: ['HG00097#1', 'HG002#1', 'HG002#2'],
    holes: [],
  })
})

// A walk file whose directory's first line is `walks`, with its header
// unit and data as walkFileOf writes them unless `data` is given
function directoryOf(walks: string, data = new Uint8Array()) {
  return WalkReader.fromFilehandles({
    data: new BlobFile(new Blob([data])),
    directory: new BlobFile(new Blob([`${walks}\n#header\t0\t0\n`])),
  })
}

test.each([
  { walks: '#walks', format: 1 },
  ...[2, 3, 4, 5, 6].map(format => ({
    walks: `#walks\tformat:i:${format}\tchunk:i:1000`,
    format,
  })),
])(
  'a walk file whose directory states $walks is refused as format $format, with what to rebuild it with',
  async ({ walks, format }) => {
    const walkFile = directoryOf(walks)
    await expect(walkFile.header()).rejects.toThrow(WalkFormatError)
    await expect(
      walkFile.cut('GRCh38#0#chr1', { start: 0, end: 600 }),
    ).rejects.toThrow(
      `This walk file is format ${format}; @jbrowse/graph-walks reads format 7, so rebuild it with gfa-to-tabix 4.0 or later`,
    )
  },
)

test.each([8, 9])(
  'a format %i walk file is refused, asking for a newer reader',
  async format => {
    await expect(
      directoryOf(`#walks\tformat:i:${format}\ttile:i:1000`).header(),
    ).rejects.toThrow(
      `This walk file is format ${format}; this version of @jbrowse/graph-walks reads format 7, so update the viewer, or @jbrowse/graph-walks, to one that reads format ${format}`,
    )
  },
)

test('a file that is not a directory is refused', async () => {
  await expect(directoryOf('GRCh38#0#chr1\t0\t100').header()).rejects.toThrow(
    "This is not a walk file's directory",
  )
})

test('a format 6 file, which has no directory, is refused by its name', async () => {
  const path = `${DATA}hprc.GRCh38.walks.bed.gz`
  await expect(
    WalkReader.open(path, url => new LocalFile(url)).header(),
  ).rejects.toThrow(
    `${path} names a walk file of format 6 or earlier, which has no .idx directory; @jbrowse/graph-walks reads format 7, so rebuild it with gfa-to-tabix 4.0 or later`,
  )
})

test('open finds the directory beside the data, before any query', async () => {
  const opened: string[] = []
  const walks = WalkReader.open(
    `${DATA}header.GRCh38.walks.gz?v=1#top`,
    url => {
      opened.push(url)
      return new LocalFile(url.replace(/[?#].*$/, ''))
    },
  )
  expect(opened).toEqual([
    `${DATA}header.GRCh38.walks.gz?v=1#top`,
    `${DATA}header.GRCh38.walks.gz.idx?v=1#top`,
  ])
  expect((await walks.header()).tile).toBe(1000)
})

// header.GRCh38's rows, with its directory changed by `edit`
async function edited(edit: (directory: string) => string) {
  const path = `${DATA}header.GRCh38.walks.gz`
  const directory = await new LocalFile(`${path}.idx`).readFile('utf8')
  return WalkReader.fromFilehandles({
    data: new LocalFile(path),
    directory: new BlobFile(new Blob([edit(directory)])),
  })
}

test('a directory and data that disagree are refused', async () => {
  await expect(
    (await edited(d => d.replace(/\tbytes:i:\d+/, '\tbytes:i:12'))).header(),
  ).rejects.toThrow(
    /directory does not match its data: the directory states 12 bytes, the data holds \d+/,
  )
  await expect(
    (await edited(d => d.replace('stay:i:311', 'stay:i:0'))).header(),
  ).rejects.toThrow(/directory does not match its data/)
  expect((await (await edited(d => d)).header()).tile).toBe(1000)
})

test('a cut keeps the reference and the haplotypes asked for', async () => {
  const names = walkNames(
    await gfa(reader('header.GRCh38'), 1000, 2600, ['HG002#2']),
  )
  expect(names).toEqual(['GRCh38#0#chr1', 'HG002#2#chr1'])
})

test('a cut lists each departure over the window, with the own steps it holds', async () => {
  const { fragments, departures } = await reader('header.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 1000, end: 2600 },
    { context: 0 },
  )
  expect(fragments.map(fragment)).toEqual([
    ['GRCh38#0#chr1', 600, 3000, '>2>3>4>5', 600, 0],
    ['HG00097#1#chr1', 600, 3000, '>2>3>4>5', 600, 0],
    ['HG002#1#chr1', 600, 2950, '>2>7>4>5', 600, 0],
    ['HG002#2#chr1', 600, 2400, '>2>3>4', 600, 20],
  ])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 1200, 550, 0, '1200+', '1800+', 'forward', 'record', '>7'],
    ['HG002#2#chr1', 2400, 20, 0, '2400+', '3000+', 'forward', 'record', '>8'],
  ])
})

// An unchopped graph whose node 1 spans 0-2600, in tile 0 alone:
//   W GRCh38 0 chr1 >1>2>4
//   W HG002  1 chr1 >1>3>4   3 replaces 2
test('a node longer than a tile is read from the tile whose extent holds it', async () => {
  const reads = await unitsRead('maxnode.GRCh38', 'GRCh38#0#chr1', walks =>
    walks.cut('GRCh38#0#chr1', { start: 2100, end: 2500 }, { context: 0 }),
  )
  expect(reads.flat()).toContain('tile 0')
  const text = await gfa(reader('maxnode.GRCh38'), 2100, 2500)
  expect(text).toMatch(/^S\t1\t/m)
  expect(text).toMatch(/^W\tGRCh38\t0\tchr1\t0\t3800\t>1>2>4$/m)
  expect(walkNames(text)).toEqual(['GRCh38#0#chr1', 'HG002#1#chr1'])
})

// Tile 0's extent is 0-1200, as node 2 runs past it; tile 1000's 1200-2400
test('a cut reads the tiles and departure blocks its window overlaps, in one read where they touch, and fetches the tile holding a leave base', async () => {
  let departures: WalkDeparture[] = []
  const reads = await unitsRead(
    'header.GRCh38',
    'GRCh38#0#chr1',
    async walks => {
      departures = (
        await walks.cut(
          'GRCh38#0#chr1',
          { start: 2100, end: 2600 },
          { context: 0 },
        )
      ).departures
    },
  )
  expect(reads).toEqual([
    ['tile 1000'],
    ['tile 2000', 'departures 1000'],
    ['tile 0'],
  ])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 1200, 550, 0, '1200+', '1800+', 'forward', 'record', '>7'],
    ['HG002#2#chr1', 2400, 20, 0, '2400+', '3000+', 'forward', 'record', '>8'],
  ])
})

// A reference that passes nodes twice, in 1 kb chunks:
//   W GRCh38 0 chr1 >1>2>3>4>5>6>5>6>7>8>7>8>9>10
//   W HG002  1 chr1 >1>2>3>4>5>6>5>6>7>8>7>11>9>10   11 replaces the second 8
// 5 and 6 are 250 bp, both passes in chunk 2000; the rest are 500 bp, so 7
// and 8 pass in chunk 3000 and again in chunk 4000
test.each([
  [2500, 3000, '5@2500 6@2750', '>5>6'],
  [4000, 5000, '7@4000 8@4500', '>7>8'],
])(
  'a node the reference passes twice is placed at its pass over the window %i-%i',
  async (start, end, placed, steps) => {
    const { graph, kept, fragments, departures } = await reader(
      'revisit.GRCh38',
    ).cut('GRCh38#0#chr1', { start, end }, { context: 0 })
    const onWindow = [...graph.nodes]
      .filter(([, node]) => node.onWindow)
      .map(([id, node]) => `${id}@${node.start}`)
    expect(onWindow.join(' ')).toBe(placed)
    expect(graphTablesGFA(graph.tables(kept, fragments, departures))).toMatch(
      new RegExp(`^W\\tGRCh38\\t0\\tchr1\\t\\d+\\t\\d+\\t${steps}$`, 'm'),
    )
  },
)

// graph-walks 4's two-argument call must not compile, or a caller that
// missed the bump would draw no departure at all
test('tables takes the departures, so the call without them is a type error', async () => {
  const { graph, kept, fragments, departures } = await reader(
    'skip.GRCh38',
  ).cut('GRCh38#0#chr1', { start: 31, end: 60 }, { context: 0 })
  // @ts-expect-error departures are required
  expect(() => graph.tables(kept, fragments)).toThrow()
  expect(graph.tables(kept, fragments, departures).departures.names).toEqual([
    'HG002#1#ctg',
  ])
})

test('a step on a second pass over the window departs from that pass', async () => {
  const { fragments, departures } = await reader('revisit.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 4000, end: 5000 },
    { context: 0 },
  )
  expect(fragments.map(fragment)).toEqual([
    ['GRCh38#0#chr1', 4000, 5000, '>7>8', 0, 0],
    ['HG002#1#chr1', 4000, 4500, '>7', 0, 500],
  ])
  expect(departures.map(departure)).toEqual([
    [
      'HG002#1#chr1',
      4500,
      500,
      0,
      '4500+',
      '5000+',
      'forward',
      'record',
      '>11',
    ],
  ])
})

// A reference that passes 7 and 8, 250 bp each, in the chunk before the
// window and again on it, with 9 and 10 between; the rest are 500 bp:
//   W GRCh38 0 chr1 >1>2>3>4>5>6>7>8>9>10>7>8>11>12
//   W HG002  1 chr1 >1>2>3>4>5>6>7>8>9>10>7>13>11>12   13 replaces the second 8
test("a cut leaves out the reference's other pass over the window's nodes", async () => {
  const { graph, kept, fragments, departures } = await reader(
    'farpass.GRCh38',
  ).cut('GRCh38#0#chr1', { start: 4500, end: 5000 }, { context: 0 })
  const text = graphTablesGFA(graph.tables(kept, fragments, departures))
  expect(text).toMatch(/^W\tGRCh38\t0\tchr1\t4500\t5000\t>7>8$/m)
  expect(text).toMatch(/^W\tHG002\t1\tchr1\t4500\t4750\t>7$/m)
  expect(text).not.toMatch(/^S\t(9|10)\t/m)
  const both = await reader('farpass.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 3900, end: 5000 },
    { context: 0 },
  )
  expect(
    graphTablesGFA(
      both.graph.tables(both.kept, both.fragments, both.departures),
    ),
  ).toMatch(/^S\t7\t\*\tLN:i:250\tSN:Z:GRCh38#0#chr1\tSO:i:4500\tSR:i:0$/m)
  expect(departures.map(departure)).toEqual([
    [
      'HG002#1#chr1',
      4750,
      250,
      0,
      '4750+',
      '5000+',
      'forward',
      'record',
      '>13',
    ],
  ])
})

// Twelve 500 bp reference nodes in 1 kb chunks (follow.GRCh38.gfa): HG002#1
// leaves chunk 1000 for 20 and 21, 2 kb each, and rejoins at 9, in chunk
// 4000; HG003#1's contig ends on 22
test('a walk leaving the window is cut there, and its record says where it rejoins', async () => {
  let cut: Awaited<ReturnType<WalkReader['cut']>> | undefined
  const reads = await unitsRead(
    'follow.GRCh38',
    'GRCh38#0#chr1',
    async walks => {
      cut = await walks.cut(
        'GRCh38#0#chr1',
        { start: 1000, end: 1500 },
        { context: 0 },
      )
    },
  )
  expect(reads).toEqual([['tile 1000', 'departures 1000']])
  const { fragments, departures } = cut!
  expect(fragments.map(fragment)).toEqual([
    ['GRCh38#0#chr1', 1000, 1500, '>3', 0, 500],
    ['HG002#1#chr1', 1000, 1500, '>3', 0, 4000],
    ['HG002#2#chr1', 1000, 1500, '>3', 0, 500],
    ['HG003#1#chr1', 1000, 1500, '>3', 0, 300],
  ])
  expect(departures.map(departure)).toEqual([
    [
      'HG002#1#chr1',
      1500,
      4000,
      0,
      '1500+',
      '4000+',
      'forward',
      'record',
      '>20>21',
    ],
  ])
})

test('the same graph in 10 kb chunks computes the departure from the steps of one piece', async () => {
  const { departures } = await reader('rejoin.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 4000, end: 4200 },
    { context: 0 },
  )
  expect(departures.map(departure)).toEqual([
    [
      'HG002#1#chr1',
      1500,
      4000,
      0,
      '1500+',
      '4000+',
      'forward',
      'steps',
      '>20>21',
    ],
  ])
})

// Twenty 500 bp reference nodes in 1 kb chunks (records.GRCh38.gfa):
//   W GRCh38 0 chr1 0    10000 >1>2>...>20
//   W HG002  1 chr1 0    2000  >1>2>3>4
//   W HG002  1 chr1 2500 3200  >30          visits no reference node
//   W HG002  1 chr1 4000 7000  >7>8>...>12
//   W HG003  1 chr1 0    10000 >1>2>3>19>20>4>5>...>18
const records = async (start: number, end: number, name: string) => {
  const { fragments, departures } = await reader('records.GRCh38').cut(
    'GRCh38#0#chr1',
    { start, end },
    { context: 0 },
  )
  return {
    fragments: fragments.filter(f => f.name === name).map(fragment),
    departures: departures.filter(d => d.name === name).map(departure),
  }
}

test("a departure across a contig's W lines counts the bp the graph lacks as missing", async () => {
  expect(await records(1000, 2000, 'HG002#1#chr1')).toEqual({
    fragments: [['HG002#1#chr1', 1000, 2000, '>3>4', 0, 0]],
    departures: [
      [
        'HG002#1#chr1',
        2000,
        2000,
        1300,
        '2000+',
        '3000+',
        'forward',
        'record',
        undefined,
      ],
    ],
  })
  expect(await records(3000, 3600, 'HG002#1#chr1')).toEqual({
    fragments: [['HG002#1#chr1', 4000, 5000, '>7>8', 0, 0]],
    departures: [
      [
        'HG002#1#chr1',
        2000,
        2000,
        1300,
        '2000+',
        '3000+',
        'forward',
        'record',
        undefined,
      ],
    ],
  })
})

test('a walk jumping away for eight chunks and back has a record each way', async () => {
  expect(await records(1000, 2000, 'HG003#1#chr1')).toEqual({
    fragments: [
      ['HG003#1#chr1', 1000, 1500, '>3', 0, 0],
      ['HG003#1#chr1', 2500, 3000, '>4', 0, 0],
    ],
    departures: [
      [
        'HG003#1#chr1',
        1500,
        0,
        0,
        '1500+',
        '9000+',
        'forward',
        'record',
        undefined,
      ],
      [
        'HG003#1#chr1',
        2500,
        0,
        0,
        '10000+',
        '1500+',
        'back',
        'record',
        undefined,
      ],
    ],
  })
})

test('a cut on a signal already aborted reads no row', async () => {
  const walks = reader('header.GRCh38')
  await walks.header()
  await walks.units('GRCh38#0#chr1')
  const read = vi.spyOn(LocalFile.prototype, 'read')
  const aborted = new AbortController()
  aborted.abort(new Error('panned away'))
  await expect(
    walks.cut(
      'GRCh38#0#chr1',
      { start: 1000, end: 2600 },
      { signal: aborted.signal },
    ),
  ).rejects.toThrow('panned away')
  expect(read).not.toHaveBeenCalled()
  read.mockRestore()
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

test('a cut over its departure budget fails without a zoom-in notice', async () => {
  const error = await reader('span.GRCh38')
    .cut('GRCh38#0#chr1', { start: 50, end: 55 }, { departureByteBudget: 1 })
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(DepartureReadError)
  expect(error).not.toBeInstanceOf(RegionTooLargeError)
  expect((error as DepartureReadError).message).toMatch(
    /^Too many departures from the reference over this window to read \([\d.]+ MB against a budget of 0\.0 MB\)/,
  )
  expect((error as DepartureReadError).limit).toBe(1)
})

// The chrM graph of HPRC's release 2, from `vg convert -f` of its GBZ and
// gfa-to-tabix --walks --refs GRCh38 --sequences
test("a cut keeps each node's bases when asked", async () => {
  const walks = reader('hprc-chrM.GRCh38')
  const { graph, kept, fragments, departures } = await walks.cut(
    'GRCh38#0#chrM',
    { start: 8000, end: 8300 },
    { sequences: true },
  )
  expect(
    walkNames(graphTablesGFA(graph.tables(kept, fragments, departures))),
  ).toHaveLength(44)
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
  expect(cut.tables.departures.names).toEqual(['HG002#2#chr1'])
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

// A walk file from rows as walkFileOf takes them, its units written in order
// or reversed, which notes the nominal start of each tile it serves
function inline(lines: string[], reversed = false, served: number[] = []) {
  const { data, directory, units } = walkFileOf(lines, reversed)
  const read = data.read.bind(data)
  data.read = (length, position = 0, opts) => {
    for (const unit of units) {
      if (
        unit.kind === 'tile' &&
        unit.offset >= position &&
        unit.offset < position + length
      ) {
        served.push(+unit.name)
      }
    }
    return read(length, position, opts)
  }
  return WalkReader.fromFilehandles({ data, directory })
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
  '#walks format:i:7 tile:i:16',
  '#reference R',
  '#haplotype H1#1',
  '#haplotype H2#1',
]
const pieces = [
  ...tinyHeader,
  'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 .',
  'R#0#chr 0 walk H1#1#c 0 0 0 2 2 .',
  'R#0#chr 0 walk H1#1#c 0 12 2 2 6 .',
  'R#0#chr 0 walk H2#1#c 0 8 1 1 4 .',
  'R#0#chr 16 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 .',
  'R#0#chr 16 walk H1#1#c 0 8 1 1 12 .',
  'R#0#chr 16 walk H1#1#c 0 20 3 4 10 .',
  'R#0#chr 16 walk H2#1#c 0 0 0 2 10 .',
  'R#0#chr 16 walk H2#1#c 0 12 2 2 14 .',
]
const stays = [
  ...tinyHeader,
  'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 .',
  'R#0#chr 0 nodes 6 4',
  'R#0#chr 0 route 4 8 -6',
  'R#0#chr 0 walk H1#1#c 0 0 0 5 2 0',
  'R#0#chr 16 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 .',
  'R#0#chr 16 nodes 2 4',
  'R#0#chr 16 route 12 -8 10',
  'R#0#chr 16 walk H1#1#c 0 20 1 4 10 .',
  'R#0#chr 16 walk H2#1#c 0 0 0 5 10 0',
]

// Each window reads both chunks, the first as the chunk before it. In
// `stays` a route's step on a reference node is the own step of a departure
// inside one piece, so the cut keeps the node, at its reference position
describe.each([
  { start: 4, end: 20, haplotypes: ['H1'], node: 6, at: 20, step: '>6' },
  { start: 12, end: 32, haplotypes: ['H2'], node: 2, at: 4, step: '>2' },
  { start: 0, end: 32, haplotypes: undefined, node: 6, at: 20, step: '>6' },
])(
  'a reference node a route steps on, cut $start-$end for $haplotypes',
  ({ start, end, haplotypes, node, at, step }) => {
    test('keeps its reference position, in whichever order the chunks come', async () => {
      const cut = (reversed?: boolean) =>
        inline(stays, reversed)
          .cutRegion({ refName: 'chr', start, end }, { haplotypes, context: 0 })
          .then(({ tables }) => ({
            text: graphTablesGFA(tables),
            own: [...tables.departures.steps].map(i => tables.nodes.names[i]),
          }))
      const { text, own } = await cut()
      expect(text).toMatch(
        new RegExp(
          `^S\\t${node}\\t\\*\\tLN:i:4\\tSN:Z:R#0#chr\\tSO:i:${at}\\tSR:i:0$`,
          'm',
        ),
      )
      expect(own).toContain(String(node))
      expect(step).toBe(`>${node}`)
      expect(await cut(true)).toEqual({ text, own })
    })
  },
)

// The pieces join into one run each, so an uncontained cut keeps each walk
// from its first step on the window to its last, the visit to the other
// chunk included, and the jumps between the pieces come from the records
test('the same walks in pieces of their own cut whole, with the jumps as records', async () => {
  const records = [
    ...pieces,
    'R#0#chr 4 24 departure H1#1#c 8 0 0 8 + 1 20 + 1',
    'R#0#chr 8 24 departure H1#1#c 12 0 0 24 + 1 8 + 1',
    'R#0#chr 4 24 departure H2#1#c 8 0 0 24 + 1 4 + 1',
    'R#0#chr 4 28 departure H2#1#c 12 0 0 8 + 1 24 + 1',
  ]
  const { tables } = await inline(records).cutRegion(
    { refName: 'chr', start: 4, end: 20 },
    { context: 0 },
  )
  const text = graphTablesGFA(tables)
  expect(text.split('\n').filter(l => l.startsWith('W'))).toEqual([
    'W\tR\t0\tchr\t4\t20\t>2>3>4>5',
    'W\tH1\t1\tc\t4\t24\t>2>6>3>4>5',
    'W\tH2\t1\tc\t0\t12\t>5>6>2',
  ])
  expect(
    [...tables.departures.leavePos].map(
      (pos, k) =>
        `${tables.departures.names[k]} ${pos}->${tables.departures.rejoinPos[k]}`,
    ),
  ).toEqual(['H1#1#c 8->20', 'H1#1#c 24->8', 'H2#1#c 24->4', 'H2#1#c 8->24'])
})

test("a route's own steps inside one piece are a departure computed from the steps", async () => {
  const { departures } = await inline(stays).cut(
    'R#0#chr',
    { start: 4, end: 20 },
    { context: 0 },
  )
  expect(departures.map(departure)).toEqual([
    ['H1#1#c', 8, 4, 0, '8+', '8+', 'adjacent', 'steps', '>6'],
  ])
})

// Four 4 bp reference nodes, the second of which the reference passes twice,
// its other pass elsewhere; H#1 runs 9 between 2 and 3:
//   R  >1>2>3>4, passes 1,2,1,1
//   H  >1>2>9>3>4
test("a departure computed from steps takes its ends' pass counts from the ref rows", async () => {
  const { graph, kept, fragments, departures } = await inline([
    '#walks format:i:7 tile:i:16 revisited:i:1',
    '#reference R',
    '#haplotype H#1',
    'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 1,2,1,1',
    'R#0#chr 0 nodes 9 4',
    'R#0#chr 0 route 4 14 -12',
    'R#0#chr 0 walk H#1#h 0 0 0 5 2 0',
  ]).cut('R#0#chr', { start: 0, end: 16 }, { context: 0 })
  expect(
    departures.map(d => [
      d.leavePos,
      d.leavePasses,
      d.rejoinPos,
      d.rejoinPasses,
      d.fromSteps,
    ]),
  ).toEqual([[8, 2, 8, 1, true]])
  expect(departures.map(placementChoice)).toEqual([true])
  const tables = graph.tables(kept, fragments, departures)
  expect(
    tables.nodes.names.map((name, i) => `${name}:${tables.nodes.passes![i]}`),
  ).toEqual(['9:0', '1:1', '2:2', '3:1', '4:1'])
  expect([...tables.departures.leavePasses]).toEqual([2])
})

test('a cut keeps to WALK_BUDGETS unless its caller lifts them', async () => {
  const window = { start: 16, end: 32 }
  const defaults = { ...WALK_BUDGETS }
  try {
    WALK_BUDGETS.departureByteBudget = 1
    await expect(
      inline(copyJump).cut('R#0#chr', window),
    ).rejects.toBeInstanceOf(DepartureReadError)
    const lifted = await inline(copyJump).cut('R#0#chr', window, {
      departureByteBudget: Infinity,
    })
    expect(lifted.departures).toHaveLength(2)
    WALK_BUDGETS.departureByteBudget = defaults.departureByteBudget
    WALK_BUDGETS.byteBudget = 1
    await expect(
      inline(copyJump).cut('R#0#chr', window),
    ).rejects.toBeInstanceOf(RegionTooLargeError)
  } finally {
    Object.assign(WALK_BUDGETS, defaults)
  }
})

test('a ref row whose passes do not count its steps is refused', async () => {
  const error = await inline([
    '#walks format:i:7 tile:i:16',
    '#reference R',
    'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4',
  ])
    .cut('R#0#chr', { start: 0, end: 16 }, { context: 0 })
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(WalkRowError)
})

// Twelve 4 bp reference nodes in 16 bp chunks; H#1 goes from 7 to 11, in
// chunk 32, and back to 6:
//   W R  0 chr >1>2>3>4>5>6>7>8>9>10>11>12
//   W H  1 h   >5>6>7>11>6>7>8
const copyJump = [
  '#walks format:i:7 tile:i:16',
  '#reference R',
  '#haplotype H#1',
  'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 .',
  'R#0#chr 16 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 .',
  'R#0#chr 16 walk H#1#h 0 0 0 3 10 .',
  'R#0#chr 16 walk H#1#h 0 16 2 3 12 .',
  'R#0#chr 32 ref R#0#chr 0 32 2 4 18,2,2,2 4,4,4,4 .',
  'R#0#chr 32 walk H#1#h 0 12 1 1 22 .',
  'R#0#chr 24 44 departure H#1#h 12 0 0 28 + 1 40 + 1',
  'R#0#chr 20 44 departure H#1#h 16 0 0 44 + 1 20 + 1',
]
const copyJumpCut = async (
  start: number,
  end: number,
  contained: boolean,
  served: number[] = [],
) => {
  const { graph, kept, fragments, departures } = await inline(
    copyJump,
    false,
    served,
  ).cut('R#0#chr', { start, end }, { context: 0, contained, nodeLimit: 5 })
  return {
    text: graphTablesGFA(graph.tables(kept, fragments, departures)),
    walks: fragments.filter(f => f.name === 'H#1#h').map(fragment),
    departures: departures.map(departure),
  }
}

test('a contained cut splits a walk where it steps onto the reference off the window', async () => {
  const detour = await copyJumpCut(20, 36, false)
  expect(detour.walks).toEqual([['H#1#h', 4, 28, '>6>7>11>6>7>8', 4, 0]])
  const records = [
    ['H#1#h', 12, 0, 0, '28+', '40+', 'forward', 'record', undefined],
    ['H#1#h', 16, 0, 0, '44+', '20+', 'back', 'record', undefined],
  ]
  expect(detour.departures).toEqual(records)
  expect(detour.text).toMatch(/^S\t11\t/m)
  expect(detour.text).toMatch(/^W\tR\t0\tchr\t20\t36\t>6>7>8>9$/m)
  const contained = await copyJumpCut(20, 36, true)
  expect(contained.walks).toEqual([
    ['H#1#h', 4, 12, '>6>7', 4, 16],
    ['H#1#h', 16, 28, '>6>7>8', 16, 0],
  ])
  expect(contained.departures).toEqual(records)
  expect(contained.text).not.toMatch(/^S\t11\t/m)
})

test('a cut reads the tiles its window overlaps, and nothing else', async () => {
  const served: number[] = []
  const { walks, departures } = await copyJumpCut(20, 28, true, served)
  expect(walks).toEqual([
    ['H#1#h', 4, 12, '>6>7', 4, 0],
    ['H#1#h', 16, 24, '>6>7', 0, 4],
  ])
  expect(departures).toHaveLength(2)
  expect(new Set(served)).toEqual(new Set([16]))
})

// Eight chunks of four 4 bp reference nodes
const eightChunks = [
  '#walks format:i:7 tile:i:16',
  '#reference R',
  ...Array.from(
    { length: 8 },
    (_, c) =>
      `R#0#chr ${16 * c} ref R#0#chr 0 ${16 * c} ${c} 4 ${8 * c + 2},2,2,2 4,4,4,4 .`,
  ),
]

test('a cut over its node limit stops decoding at the tile that puts it over', async () => {
  const error = await inline(eightChunks)
    .cutRegion(
      { refName: 'chr', start: 0, end: 128 },
      { context: 0, nodeLimit: 6 },
    )
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(RegionTooLargeError)
  expect(error).toMatchObject({ limit: 6, windowBp: 128, fitsBp: 24 })
  const whole = await inline(eightChunks).cutRegion(
    { refName: 'chr', start: 0, end: 128 },
    { context: 0, nodeLimit: 32 },
  )
  expect(whole.tables.nodes.names).toHaveLength(32)
})

// gfa-to-tabix's walks-two-refs.gfa on R in 10 bp chunks, once as written and
// once with --abridge-bytes 1, which abridges every walk in every chunk
async function twoRefs(name: string, expand?: (name: string) => boolean) {
  const { graph, kept, fragments, departures } = await reader(name).cut(
    'R#0#chr',
    { start: 0, end: 30 },
    { expand },
  )
  return {
    tables: graph.tables(kept, fragments, departures),
    fragments,
    departures,
  }
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
  const abridged = new Set(fragments.filter(f => f.abridged).map(f => f.name))
  expect([...abridged].sort()).toEqual(['D#1#d', 'H#1#h', 'I#1#i', 'Q#0#chr'])
  names.forEach((name, i) => {
    expect(Number.isNaN(abridgedSteps![i])).toBe(
      name === 'R#0#chr' ||
        !abridged.has(name) ||
        fragments.some(
          f =>
            f.name === name &&
            !f.abridged &&
            f.hapStart === tables.walks.starts[i],
        ),
    )
  })
})

test('an expanded walk is read from its detail records as the unabridged file cuts it', async () => {
  const all = () => true
  const expanded = await twoRefs('abridged.R', all)
  const whole = await twoRefs('unabridged.R')
  expect(graphTablesGFA(expanded.tables)).toBe(graphTablesGFA(whole.tables))
  expect(expanded.departures.map(departure)).toEqual(
    whole.departures.map(departure),
  )
  expect(whole.departures.map(departure)).toEqual([
    ['D#1#d', 24, 0, 0, '24+', '0+', 'back', 'record', undefined],
    ['H#1#h', 118, 3, 0, '18+', '18+', 'adjacent', 'steps', '>7'],
    ['I#1#i', 6, 0, 0, '6+', '18-', 'switch', 'record', undefined],
    ['I#1#i', 18, 0, 0, '6-', '18+', 'switch', 'record', undefined],
    ['Q#0#chr', 6, 6, 0, '6+', '12+', 'forward', 'record', '>6'],
  ])
})

test('the reference sequences leave out detail records and departures', async () => {
  expect(await reader('abridged.R').referenceSequences()).toEqual(['R#0#chr'])
})

// abridged.R, each read of a detail unit made through `around`
async function heldDetails(
  around: <T>(
    read: () => Promise<T>,
    detail: { start: number; signal?: AbortSignal },
  ) => Promise<T>,
) {
  const path = `${DATA}abridged.R.walks.gz`
  const units = await reader('abridged.R').units('R#0#chr')
  const details = new Set(
    units.filter(u => u.kind === 'detail').map(u => u.offset),
  )
  const data = new LocalFile(path)
  const read = data.read.bind(data)
  data.read = (length, position = 0, opts = {}) =>
    details.has(position)
      ? around(() => read(length, position, opts), {
          start: position,
          signal: opts.signal,
        })
      : read(length, position, opts)
  return WalkReader.fromFilehandles({
    data,
    directory: new LocalFile(`${path}.idx`),
  })
}

async function expandedCut(walks: WalkReader, signal?: AbortSignal) {
  const { graph, kept, fragments, departures } = await walks.cut(
    'R#0#chr',
    { start: 0, end: 30 },
    { expand: () => true, signal },
  )
  return {
    gfa: graphTablesGFA(graph.tables(kept, fragments, departures)),
    fragments,
  }
}

test('detail records read at once, finishing in reverse, cut as reads one at a time do', async () => {
  const pieces = (await twoRefs('abridged.R')).fragments.flatMap(f =>
    f.abridged
      ? [`${f.abridged.detail.refName} ${f.abridged.detail.name}`]
      : [],
  )
  const details = new Set(pieces).size
  expect(details).toBeLessThan(pieces.length)
  let queue = Promise.resolve()
  let reading = 0
  let most = 0
  const serial = await expandedCut(
    await heldDetails(read => {
      const turn = queue.then(async () => {
        most = Math.max(most, ++reading)
        const bytes = await read()
        reading--
        return bytes
      })
      queue = turn.then(() => undefined)
      return turn
    }),
  )
  expect(most).toBe(1)
  const waiting: (() => void)[] = []
  const finished: number[] = []
  const concurrent = await expandedCut(
    await heldDetails(async (read, { start }) => {
      const bytes = await read()
      await new Promise<void>(resolve => {
        waiting.push(() => {
          finished.push(start)
          resolve()
        })
        if (waiting.length === details) {
          for (const release of waiting.toReversed()) {
            release()
          }
        }
      })
      return bytes
    }),
  )
  expect(finished).toHaveLength(details)
  expect(finished).not.toEqual(finished.toSorted((a, b) => a - b))
  expect(concurrent).toEqual(serial)
})

test('aborting an expanding cut stops its detail reads', async () => {
  const controller = new AbortController()
  const signals: AbortSignal[] = []
  const walks = await heldDetails(
    (_read, { signal }) =>
      new Promise<never>((_resolve, reject) => {
        signals.push(signal!)
        signal!.addEventListener('abort', () => {
          reject(signal!.reason as Error)
        })
        if (signals.length === 3) {
          controller.abort(new Error('panned away'))
        }
      }),
  )
  await expect(expandedCut(walks, controller.signal)).rejects.toThrow(
    'panned away',
  )
  expect(signals).toHaveLength(3)
  expect(signals.every(signal => signal.aborted)).toBe(true)
})

test('readInOrder keeps at most its limit reading, and uses results in order', async () => {
  let reading = 0
  let most = 0
  const used: number[] = []
  await readInOrder(
    [5, 1, 4, 2, 3, 0],
    3,
    new AbortController().signal,
    async (ms: number) => {
      most = Math.max(most, ++reading)
      await new Promise(resolve => setTimeout(resolve, ms))
      reading--
      return ms * 10
    },
    (_ms, result) => used.push(result),
  )
  expect(most).toBe(3)
  expect(used).toEqual([50, 10, 40, 20, 30, 0])
})

// gfa-to-tabix --walks --chunk 10 of >1>2>9>10 against a reference >1 to >10
// of 10 bp nodes: HG002#1#ctg leaves the reference after 20 and rejoins at 80,
// and its record is filed over the stretch between
// The drawn stretch is the window and its context, so a cut keeps the same
// reference whether the context is asked for or the window widened by it, on
// both sides and past the chunks under the window itself
test('a cut reads the chunks under its context, not only under its window', async () => {
  const walks = reader('skip.GRCh38')
  const reference = async (start: number, end: number, context: number) =>
    (
      await walks.cut('GRCh38#0#chr1', { start, end }, { context })
    ).fragments.find(f => f.name === 'GRCh38#0#chr1')?.refStarts
  expect(await reference(32, 38, 5)).toEqual(await reference(27, 43, 0))
  expect(await reference(32, 38, 15)).toEqual(await reference(17, 53, 0))
})

test('a window inside a jump finds the jump in the departure read', async () => {
  const walks = reader('skip.GRCh38')
  const { fragments, departures } = await walks.cut(
    'GRCh38#0#chr1',
    { start: 31, end: 60 },
    { context: 0 },
  )
  expect(fragments.map(f => f.name)).toEqual(['GRCh38#0#chr1'])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#ctg', 20, 0, 0, '20+', '80+', 'forward', 'record', undefined],
  ])
  const { tables } = await walks.cutRegion({
    refName: 'GRCh38#0#chr1',
    start: 31,
    end: 60,
  })
  expect(tables.departures.names).toEqual(['HG002#1#ctg'])
  expect([
    ...tables.departures.leavePos,
    ...tables.departures.rejoinPos,
  ]).toEqual([20, 80])
  const filtered = await walks.cut(
    'GRCh38#0#chr1',
    { start: 31, end: 60 },
    { keep: name => name === 'GRCh38#0#chr1' },
  )
  expect(filtered.departures).toEqual([])
})

// gfa-to-tabix --walks --chunk 20 of test_data/span.gfa, ten 10 bp reference
// nodes and three off it:
//   W GRCh38 0 chr1 >1>2>3>4>5>6>7>8>9>10
//   W HG002  1 chr1 >1>2>3>4>20>21>22>8>9>10   its own steps filed in chunk 20
//   W HG003  1 chr1 <10<9<8<22<21<20<4<3<2<1   its own steps filed in chunk 60
//   W HG004  1 chr1 <10<9<8<5<4<3<2<1          jumps from 70 to 50
// A cut of 50-55 reads chunks 20 and 40
const spanCut = (haplotype: string, opts = {}, walks = reader('span.GRCh38')) =>
  walks.cut(
    'GRCh38#0#chr1',
    { start: 50, end: 55 },
    { context: 0, keep: walkNameFilter([haplotype], 'GRCh38'), ...opts },
  )

test("a record under the attach threshold carries the walk's own steps from the read", async () => {
  const { fragments, departures } = await spanCut('HG002#1')
  expect(fragments.map(f => f.name)).toEqual(['GRCh38#0#chr1'])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 40, 30, 0, '40+', '70+', 'forward', 'record', '>20>21>22'],
  ])
})

test('a leave piece the read lacks is fetched for the own steps, once per tile, within the fetch cap', async () => {
  let departures: WalkDeparture[] = []
  const reads = await unitsRead('span.GRCh38', 'GRCh38#0#chr1', async walks => {
    departures = (await spanCut('HG003#1', {}, walks)).departures
  })
  expect(departures.map(departure)).toEqual([
    ['HG003#1#chr1', 30, 30, 0, '70-', '40-', 'forward', 'record', '<22<21<20'],
  ])
  expect(reads).toEqual([['tile 40', 'departures 20'], ['tile 60']])
  const capped = await spanCut('HG003#1', { attachFetches: 0 })
  expect(capped.departures.map(departure)).toEqual([
    ['HG003#1#chr1', 30, 30, 0, '70-', '40-', 'forward', 'record', undefined],
  ])
})

test('a record over the attach threshold carries no steps', async () => {
  const { departures, kept } = await spanCut('HG002#1', { attachBp: 29 })
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 40, 30, 0, '40+', '70+', 'forward', 'record', undefined],
  ])
  expect([...kept].sort((a, b) => a - b)).toEqual([6])
})

test('a direct jump across the window is a record with nothing of its own', async () => {
  const { fragments, departures } = await spanCut('HG004#1')
  expect(fragments.map(f => f.name)).toEqual(['GRCh38#0#chr1'])
  expect(departures.map(departure)).toEqual([
    ['HG004#1#chr1', 30, 0, 0, '70-', '50-', 'forward', 'record', undefined],
  ])
})

// Twelve 10 bp reference nodes in 20 bp chunks (attach.GRCh38.gfa); HG003#1
// leaves after 2, in chunk 0, for 21 and HG002#1 after 3, in chunk 20, for
// 20, and both rejoin at 11, in chunk 100:
//   W GRCh38 0 chr1 >1>2>3>4>5>6>7>8>9>10>11>12
//   W HG002  1 chr1 >1>2>3>20>11>12
//   W HG003  1 chr1 >1>2>21>11>12
test('two records whose leave pieces lie in adjacent chunks both get their own steps', async () => {
  const { departures } = await reader('attach.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 100, end: 120 },
    { context: 0 },
  )
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 30, 10, 0, '30+', '100+', 'forward', 'record', '>20'],
    ['HG003#1#chr1', 20, 10, 0, '20+', '100+', 'forward', 'record', '>21'],
  ])
})

test('a table entry of length 0 is never read', async () => {
  const { data, directory, units } = walkFileOf(copyJump)
  const bytes = new Uint8Array(await data.readFile())
  const end = bytes.length - 28
  const entries = [
    ['#unit', 'R#0#chr', 0, 48, 'table', '.'],
    ['#entry', 'R#0#chr', 0, 48, 'tile', '0', 0, 0],
    ...units.map(u => [
      '#entry',
      u.refName,
      u.start,
      u.end,
      u.kind,
      u.name,
      u.offset,
      u.length,
    ]),
  ]
  const table = bgzf(`${entries.map(e => e.join('\t')).join('\n')}\n`)
  const tableLength = table.reduce((sum, block) => sum + block.length, 0)
  const parts = [bytes.subarray(0, end), ...table, bytes.subarray(end)]
  const size = parts.reduce((sum, part) => sum + part.length, 0)
  const [walks = '', header = ''] = (await directory.readFile('utf8')).split(
    '\n',
  )
  const index = `${walks.replace(/\tbytes:i:\d+/, `\tbytes:i:${size}`)}\n${header}\nR#0#chr\t${end}\t${tableLength}\n`
  const reader = WalkReader.fromFilehandles({
    data: new BlobFile(new Blob(parts)),
    directory: new BlobFile(new Blob([index])),
  })
  expect((await reader.units('R#0#chr'))[0]).toMatchObject({ length: 0 })
  const { departures } = await reader.cut(
    'R#0#chr',
    { start: 0, end: 48 },
    { context: 0 },
  )
  expect(departures).toHaveLength(2)
})

test("a sequence's length, the bytes a cut reads and the reference nodes over a window come from its tiles", async () => {
  const walks = reader('header.GRCh38')
  const units = await walks.units('GRCh38#0#chr1')
  const tile = (name: string) =>
    units.find(u => u.kind === 'tile' && u.name === name)!
  expect(await walks.sequenceLength('GRCh38#0#chr1')).toBe(3600)
  expect(await walks.sequenceLength('GRCh38#0#chr2')).toBe(0)
  expect(
    await walks.bytesFor(
      'GRCh38#0#chr1',
      { start: 2100, end: 2600 },
      {
        context: 0,
      },
    ),
  ).toBe(
    tile('1000').length +
      tile('2000').length +
      units.find(u => u.kind === 'departures' && u.start === 2000)!.length,
  )
  expect(
    (
      await walks.referenceNodes('GRCh38#0#chr1', { start: 1100, end: 2500 })
    ).map(n => `${n.id}@${n.start}-${n.end}`),
  ).toEqual(['2@600-1200', '3@1200-1800', '4@1800-2400', '5@2400-3000'])
})
