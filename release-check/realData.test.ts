import { TabixIndexedFile } from '@gmod/tabix'
import { RemoteFile } from 'generic-filehandle2'
import { expect, test } from 'vitest'

import { formContradictions } from '../src/formsCheck.ts'
import { structuralForms } from '../src/structuralForms.ts'
import { formsCut } from '../src/walkDepartures.ts'
import { WalkReader } from '../src/walkReader.ts'

import type { HaplotypeDepartures } from '../src/walkDepartures.ts'

// The 28 MC chr22 windows of gfa-to-tabix's
// docs/departures/results/minimal/mc22.forms-compare-review.txt
const WINDOWS = [
  [18_600_000, 18_680_000],
  [21_500_000, 21_600_000],
  [23_000_000, 23_100_000],
  [30_000_000, 30_010_000],
  [42_100_000, 42_150_000],
  [50_000_000, 50_100_000],
  [18_900_000, 19_000_000],
  [19_000_000, 19_100_000],
  [20_300_000, 20_400_000],
  [20_700_000, 20_800_000],
  [21_000_000, 21_100_000],
  [21_400_000, 21_500_000],
  [18_700_000, 18_750_000],
  [22_000_000, 22_100_000],
  [22_300_000, 22_400_000],
  [22_500_000, 22_600_000],
  [22_800_000, 22_900_000],
  [22_000_000, 22_900_000],
  [12_345_000, 12_395_000],
  [15_000_000, 15_050_000],
  [25_500_000, 25_550_000],
  [33_333_000, 33_383_000],
  [38_000_000, 38_050_000],
  [45_678_000, 45_728_000],
  [48_000_000, 48_050_000],
  [10_500_000, 10_550_000],
  [50_780_000, 50_818_468],
  [30_000_000, 30_001_000],
] as const

const walks = process.env.WALKS ?? ''
const refName = process.env.REF_NAME ?? 'chr22'
const lifted = {
  stepBudget: Infinity,
  byteBudget: Infinity,
  departureByteBudget: Infinity,
}

const reader = () =>
  new WalkReader(
    walks.startsWith('http')
      ? new TabixIndexedFile({
          filehandle: new RemoteFile(walks),
          tbiFilehandle: new RemoteFile(`${walks}.tbi`),
        })
      : new TabixIndexedFile({ path: walks, tbiPath: `${walks}.tbi` }),
  )

const net = (h: HaplotypeDepartures) =>
  h.departures.reduce((s, d) => s + d.altLength - (d.refEnd - d.refStart), 0)

const summary = (h: HaplotypeDepartures) =>
  `${net(h)} ${JSON.stringify(h.covered)}`

test.runIf(walks !== '').each(WINDOWS)(
  '%i-%i: abridged pieces read unexpanded as expanded, context on both sides, and no form contradicts a member',
  async (start, end) => {
    const region = { refName, start, end }
    const cut = (expand: boolean) =>
      formsCut(reader(), region, { ...lifted, expand, partial: false })
    const [unexpanded, expanded] = [await cut(false), await cut(true)]
    expect(new Map(unexpanded.map(h => [h.haplotype, summary(h)]))).toEqual(
      new Map(expanded.map(h => [h.haplotype, summary(h)])),
    )

    const sequence =
      (await reader().sequenceFor(
        (await reader().header()).references[0] ?? '',
        refName,
      )) ?? refName
    const reference = async (from: number, to: number, context: number) =>
      (
        await reader().cut(
          sequence,
          { start: from, end: to },
          { ...lifted, context },
        )
      ).fragments.find(f => f.name === sequence)?.refStarts
    expect(await reference(start, end, 1000)).toEqual(
      await reference(start - 1000, end + 1000, 0),
    )

    const window = { start, end }
    for (const floorBp of ['auto', 50, 1000, 10_000] as const) {
      for (const foldRareForms of [true, false]) {
        const result = structuralForms(expanded, {
          window,
          floorBp,
          foldRareForms,
          exclude: ['CHM13', 'GRCh38'],
        })
        expect(formContradictions(expanded, result, { window })).toEqual([])
      }
    }
  },
)
