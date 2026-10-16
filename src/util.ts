// A cut's tables: what a GFA of S lines tagged LN/SN/SO/SR, L lines and W
// lines states, every link end and walk step an index into the node table.
// @jbrowse/bandage-core's GraphTables is the same shape.
export interface GraphTables {
  nodes: {
    names: string[]
    lengths: Int32Array
    refs: Int32Array
    starts: Float64Array
    ranks: Int32Array
    refNames: string[]
  }
  links: {
    from: Int32Array
    to: Int32Array
    strands: Uint8Array
  }
  walks: {
    names: string[]
    starts: Float64Array
    ends: Float64Array
    offsets: Int32Array
    steps: Int32Array
    reversed: Uint8Array
    // each fragment's W record start, and why it stops at either end, from a
    // walk-file cut; absent from tables another source builds
    fragStarts?: Float64Array
    before?: WalkStopTable
    after?: WalkStopTable
  }
}

// The kinds a WalkStop names, coded by index
export const WALK_STOP_KINDS = [
  'window',
  'elsewhere',
  'break',
  'end',
  'recordEnd',
] as const

// A WalkStop per walk: its kind's index in WALK_STOP_KINDS, `bp` NaN where the
// stop states none, and the chunk as an index into `chunkRefNames`, -1 for none
export interface WalkStopTable {
  kinds: Uint8Array
  bp: Float64Array
  chunkRefs: Int32Array
  chunkStarts: Float64Array
  chunkRefNames: string[]
}

/**
 * A window too large to cut. `regionTooLarge` and the name are what JBrowse
 * and @jbrowse/bandage-core read a zoom-in notice by; `fitsBp` is about how
 * wide a window would fit.
 */
export class RegionTooLargeError extends Error {
  readonly regionTooLarge = true
  override name = 'NodeLimitError'

  readonly limit: number
  readonly windowBp: number
  readonly fitsBp: number

  constructor(
    message: string,
    limit: number,
    windowBp: number,
    fitsBp: number,
  ) {
    super(message)
    this.limit = limit
    this.windowBp = windowBp
    this.fitsBp = fitsBp
  }
}

const SEP = '#'

// `grape#1#chr1` -> `grape`
export function panSNSample(refName: string) {
  return refName.split(SEP)[0]!
}

// `grape#1#chr1` and `grape#chr1` -> `chr1`
export function panSNContig(refName: string) {
  const parts = refName.split(SEP)
  return parts.length >= 3
    ? parts.slice(2).join(SEP)
    : parts.length === 2
      ? parts[1]!
      : refName
}

export function panSNHaplotype(refName: string) {
  const parts = refName.split(SEP)
  return parts.length >= 3 ? `${parts[0]}${SEP}${parts[1]}` : undefined
}

export function panSNMatchesPrefix(
  refName: string,
  prefix: string | undefined,
) {
  return (
    prefix !== undefined &&
    (refName === prefix || refName.startsWith(prefix + SEP))
  )
}

const precise = (n: number) => Number.parseFloat(n.toPrecision(3))

// `1.5Kbp`, `12Mbp`, `850bp`, as JBrowse writes a span
export function bpString(bp: number) {
  const mb = precise(bp / 1_000_000)
  return mb >= 1
    ? `${mb.toLocaleString('en-US')}Mbp`
    : bp >= 1000
      ? `${precise(bp / 1000).toLocaleString('en-US')}Kbp`
      : `${Math.floor(bp)}bp`
}
