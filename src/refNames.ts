import { panSNContig, panSNHaplotype, panSNSample } from './util.ts'

// A tab can occur in neither half: both come out of BED columns
function qualifiedKey(assemblyName: string, refName: string) {
  return `${assemblyName}\t${refName}`
}

/**
 * A file's sequence names under each spelling a caller may ask by: the file's
 * own (`GRCh38#0#chr6`), and a PanSN name's contig under its sample (`chr6`
 * on `GRCh38`) and under its haplotype (`chr6` on `GRCh38#0`). A pggb name
 * with no haplotype field, `chm13#chr1`, is `chr1` on `chm13`.
 */
export function refNameLookup(names: string[]) {
  const lookup = new Map<string, string>()
  for (const name of names) {
    lookup.set(name, name)
    const sample = panSNSample(name)
    if (sample !== name) {
      const contig = panSNContig(name)
      lookup.set(qualifiedKey(sample, contig), name)
      const haplotype = panSNHaplotype(name)
      if (haplotype !== undefined) {
        lookup.set(qualifiedKey(haplotype, contig), name)
      }
    }
  }
  return lookup
}

/**
 * The file's name for `refName` on `assemblyName`, a sample or a haplotype,
 * else for `refName` as the file spells it; undefined when it holds neither
 */
export function resolveRefName(
  lookup: Map<string, string>,
  assemblyName: string,
  refName: string,
) {
  return lookup.get(qualifiedKey(assemblyName, refName)) ?? lookup.get(refName)
}
