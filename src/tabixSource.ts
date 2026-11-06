/**
 * A Tabix-indexed file as the node and link index is read: @gmod/tabix's
 * TabixIndexedFile is one. `bytesForRegions` estimates the compressed bytes a
 * read fetches.
 */
export interface TabixSource {
  getHeaderLines(opts?: { signal?: AbortSignal }): Promise<string[]>
  getLines(
    refName: string,
    start: number,
    end: number,
    opts: { signal?: AbortSignal; lineCallback: (line: string) => void },
  ): Promise<void>
  getReferenceSequenceNames(opts?: { signal?: AbortSignal }): Promise<string[]>
  bytesForRegions?(
    regions: { refName: string; start: number; end: number }[],
    opts?: { signal?: AbortSignal },
  ): Promise<number>
}

// @gmod/tabix loops forever on a NaN end and reads nothing for a NaN start
// (GMOD/tabix-js#157), so every query refuses a range that is not finite first
export function checkRange(refName: string, start: number, end: number) {
  if (!Number.isFinite(start) || !Number.isFinite(end)) {
    throw new Error(
      `${refName}:${start}-${end} is not a finite range to query the graph's index for`,
    )
  }
}
