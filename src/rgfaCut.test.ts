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

test('a cut is byte-identical across calls', async () => {
  const cut = () => rgfaCut(files('rgfa_ecoli'), window, { hops: 1 })
  expect(await cut()).toBe(await cut())
})

test('a cut refuses a range that is not finite before querying', async () => {
  await expect(
    rgfaCut(files('rgfa_ecoli'), { ...window, end: Number.NaN }),
  ).rejects.toThrow("is not a finite range to query the graph's index for")
})
