import { TabixIndexedFile } from '@gmod/tabix'

import { carrierNames, parseGraphHeader, parseGraphRow } from './rgfaBed.ts'
import { rgfaCut } from './rgfaCut.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

const open = (path: string) =>
  new TabixIndexedFile({ path, tbiPath: `${path}.tbi` })

// gfa-to-tabix's ecoli_pggb_subgraph.gfa on K12, written in both formats
const formats = (layout: string) => ({
  one: {
    segments: open(`${DATA}pggb_ecoli_${layout}.segs.bed.gz`),
    links: open(`${DATA}pggb_ecoli_${layout}.links.bed.gz`),
  },
  two: { graph: open(`${DATA}pggb_ecoli_${layout}.graph.bed.gz`) },
})

// The cut's lines with each SM:Z: list as a sorted set, since format 2 names
// carriers in #samples order
function lines(gfa: string) {
  return gfa.split('\n').map(line =>
    line
      .split('\t')
      .map(col =>
        col.startsWith('SM:Z:')
          ? `SM:Z:${col.slice(5).split(',').sort().join(',')}`
          : col,
      )
      .join('\t'),
  )
}
const withoutSamples = (gfa: string) =>
  lines(gfa).map(l => l.replace(/\tSM:Z:\S*/, ''))

test.each(['anchored', 'contig'])(
  'a %s cut of format 2 holds format 1 cut',
  async layout => {
    const { one, two } = formats(layout)
    for (let start = 1004500; start < 1004961; start += 37) {
      const window = { refName: 'K12#1#chr', start, end: start + 60 }
      const a = await rgfaCut(one, window, { hops: 1 })
      const b = await rgfaCut(two, window, { hops: 1 })
      expect(withoutSamples(b)).toEqual(withoutSamples(a))
      // every carrier list format 2 gives matches format 1's
      const named = new Set(lines(a))
      for (const line of lines(b).filter(l => l.includes('SM:Z:'))) {
        expect(named).toContain(line)
      }
    }
  },
)

test('a format 2 node row names its carriers from #samples', () => {
  const header = parseGraphHeader([
    '#gfa-to-tabix\tformat:i:2\tlayout:Z:anchored\tgraph:Z:7b113ffb87119fe1',
    '#samples\tCFT073.1,IAI39.1,K12.1,NCTC86.1,Sakai.1',
  ])
  expect(header).toEqual({
    anchored: true,
    samples: ['CFT073.1', 'IAI39.1', 'K12.1', 'NCTC86.1', 'Sakai.1'],
  })
  expect(carrierNames('c1', header.samples)).toEqual([
    'K12.1',
    'NCTC86.1',
    'Sakai.1',
  ])
  const row = parseGraphRow(
    'K12#1#chr\t1004500\t1004501\tN\t1\t0\tK12#1#chr\t1004500\t1004501\tSB:H:c',
    header,
  )
  expect(row).toEqual({
    kind: 'N',
    segment: {
      refName: 'K12#1#chr',
      start: 1004500,
      end: 1004501,
      id: '1',
      rank: 0,
      tags: 'SM:Z:K12.1,NCTC86.1',
    },
  })
  expect(parseGraphRow('K12#1#chr\t0\t1\tX\tsomething', header)).toBe(undefined)
})

test('a format other than 2 is refused', () => {
  expect(() => parseGraphHeader(['#gfa-to-tabix\tformat:i:3'])).toThrow(
    'gfa-to-tabix format 3: this graph-walks reads format 2',
  )
  expect(() => parseGraphHeader(['#other'])).toThrow('no #gfa-to-tabix header')
})
