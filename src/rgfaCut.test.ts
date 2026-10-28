import { TabixIndexedFile } from '@gmod/tabix'

import { rgfaCut } from './rgfaCut.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

function files(name: string) {
  const open = (path: string) =>
    new TabixIndexedFile({ path, tbiPath: `${path}.tbi` })
  return {
    segments: open(`${DATA}${name}.segs.bed.gz`),
    links: open(`${DATA}${name}.links.bed.gz`),
  }
}

const segmentIds = (gfa: string) =>
  gfa
    .split('\n')
    .filter(l => l.startsWith('S\t'))
    .map(l => l.split('\t')[1]!)

// The paa window of an E. coli rGFA, in the contig layout (`rgfa_ecoli`) and
// in gfa-to-tabix's anchored layout (`rgfa_ecoli_anchored`)
const window = { refName: 'K12#1#chr', start: 993236, end: 997574 }

test('a contig-layout cut reaches off-reference neighbours as it hops', async () => {
  const none = segmentIds(await rgfaCut(files('rgfa_ecoli'), window))
  const one = segmentIds(
    await rgfaCut(files('rgfa_ecoli'), window, { hops: 1 }),
  )
  expect(one).toEqual(expect.arrayContaining(none))
  expect(one.length).toBeGreaterThan(none.length)
})

test('an anchored cut holds whole without hops what a contig cut hops to', async () => {
  const walked = segmentIds(
    await rgfaCut(files('rgfa_ecoli'), window, { hops: 20 }),
  )
  const anchored = await rgfaCut(files('rgfa_ecoli_anchored'), window)
  expect(segmentIds(anchored)).toEqual(expect.arrayContaining(walked))
  expect(await rgfaCut(files('rgfa_ecoli_anchored'), window, { hops: 3 })).toBe(
    anchored,
  )
})

// gfa-to-tabix's tests/data/rgfa.gfa indexed with `--max-component 0`, which
// leaves s3's one-node bubble on its own coordinates
test('an anchored cut marks the segments of a moved component', async () => {
  const gfa = await rgfaCut(files('rgfa_moved'), {
    refName: 'chr1',
    start: 0,
    end: 13,
  })
  expect(gfa.split('\n').filter(l => l.includes('mv:i:1'))).toEqual([
    'S\ts3\t*\tLN:i:3\tSN:Z:b#1#c\tSO:i:100\tSR:i:1\tmv:i:1',
  ])
  expect(gfa).toContain('L\ts1\t+\ts3\t-\t0M')
})

test('a cut is byte-identical across calls', async () => {
  const cut = () => rgfaCut(files('rgfa_ecoli'), window, { hops: 1 })
  expect(await cut()).toBe(await cut())
})

test('a cut refuses a range that is not finite before querying', async () => {
  await expect(
    rgfaCut(files('rgfa_ecoli'), { ...window, end: Number.NaN }),
  ).rejects.toThrow("is not a finite range to query the graph's index for")
})
