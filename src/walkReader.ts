import { refNameLookup, resolveRefName } from './refNames.ts'
import { decodeChunk, pieceLinks } from './walkChunk.ts'
import {
  WalkGraph,
  byteBudgetError,
  chunkQueryStart,
  joinPieces,
  keptNodesError,
  keptStepsError,
  lookbackChunks,
  stepBudgetError,
  walkCut,
  walkHeader,
  walksKept,
} from './walkRows.ts'

import type { DecodedChunk } from './walkChunk.ts'
import type { BudgetNames, WalkHeader, WalkRow } from './walkRows.ts'

// A window reads whole chunks and keeps the steps on it, so a narrow window
// reads several times what it keeps: 464 haplotypes across the 1.5 kb ABCA7
// repeat read 5.0 M steps and keep 1.16 M, and the cut takes 0.4 s from warm
// files. The step budget bounds what the cut keeps, and rows are read to this
// multiple of it.
const READ_STEPS_PER_KEPT = 8

// A walk the read cuts off goes on in the chunk its row names, out to which a
// cut reads on, in at most this many rounds and this many chunks from the read
const FOLLOW_ROUNDS = 3
const FOLLOW_CHUNKS = 4

// bp of reference a cut keeps past each end of the window
const WALK_CONTEXT = 1000

// The budgets the JBrowse graph plugin and BandageJS cut a browser's window by
export const WALK_BUDGETS = { stepBudget: 4_000_000, byteBudget: 8_000_000 }

/**
 * A Tabix-indexed file as the reader uses one: @gmod/tabix's
 * TabixIndexedFile is one. `bytesForRegions` estimates the compressed bytes a
 * read fetches, which the byte budget needs.
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

// A walk file in a format other than 3, which the reader refuses
export class WalkFormatError extends Error {
  override name = 'WalkFormatError'
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

// Settles like `promise` unless `signal` aborts first
function raceAbort<T>(promise: Promise<T>, signal?: AbortSignal) {
  if (!signal) {
    return promise
  }
  signal.throwIfAborted()
  return new Promise<T>((resolve, reject) => {
    const onAbort = () => {
      reject(
        signal.reason instanceof Error
          ? signal.reason
          : new Error(String(signal.reason)),
      )
    }
    signal.addEventListener('abort', onAbort, { once: true })
    promise.then(resolve, reject).finally(() => {
      signal.removeEventListener('abort', onAbort)
    })
  })
}

export interface WalkCutOptions {
  signal?: AbortSignal
  // which walks to decode, by `sample#haplotype#contig`; undefined is all
  keep?: (name: string) => boolean
  // the most steps a cut keeps; Infinity, the default, has no limit
  stepBudget?: number
  // the most compressed bytes a cut fetches; Infinity, the default, has no
  // limit
  byteBudget?: number
  context?: number
  // false lists no links, for a caller that draws no edges
  links?: boolean
  // keep the bases of the nodes, from a file written with --sequences
  sequences?: boolean
  // how a zoom-in notice names the budgets
  budgetNames?: BudgetNames
}

export interface RegionCutOptions extends Omit<WalkCutOptions, 'keep'> {
  // the `sample#haplotype` walks to keep beside the reference; every walk when
  // undefined or empty
  haplotypes?: string[]
  // the sample whose sequence the region names, the file's first reference
  // by default
  referenceSample?: string
  // the most nodes a cut keeps; Infinity, the default, has no limit
  nodeLimit?: number
}

export type WalkReaderHeader = WalkHeader & { chunk: number; lookback: number }

/**
 * A `gfa-to-tabix --walks` file, `<prefix>.<sample>.walks.bed.gz`, read a
 * window at a time: one read over whole chunks, from the lookback chunks
 * before the window to its end, each chunk decoded once. The index first
 * estimates the bytes the read fetches, and a window past the byte budget
 * reads no row. Only the walks `keep` accepts are decoded. The step budget
 * bounds the steps the cut keeps, and rows holding READ_STEPS_PER_KEPT times
 * as many stop the read.
 */
export class WalkReader {
  private headerRead: Promise<WalkReaderHeader> | undefined
  private lookupRead: Promise<Map<string, string>> | undefined
  private readonly walks: TabixSource

  constructor(walks: TabixSource) {
    this.walks = walks
  }

  // The header, read once for every caller: one caller giving up stops only
  // its own wait, and a read that fails is tried again by the next caller
  header(opts: { signal?: AbortSignal } = {}) {
    if (!this.headerRead) {
      const read = this.walks.getHeaderLines().then(lines => {
        const header = walkHeader(lines)
        if (header.format !== 3) {
          throw new WalkFormatError(
            `This walk file is format ${header.format ?? 1}; @jbrowse/graph-walks reads format 3, so rebuild it with gfa-to-tabix 1.1 or later`,
          )
        }
        const chunk = header.chunk ?? 65536
        return {
          ...header,
          chunk,
          lookback: lookbackChunks(header.maxNode, chunk),
        }
      })
      read.catch(() => {
        if (this.headerRead === read) {
          this.headerRead = undefined
        }
      })
      this.headerRead = read
    }
    return raceAbort(this.headerRead, opts.signal)
  }

  // The reference sequences the walk file holds windows on
  referenceSequences(opts: { signal?: AbortSignal } = {}) {
    return this.walks.getReferenceSequenceNames(opts)
  }

  /**
   * The walk file's sequence for `refName` on `sample`, or for `refName` as
   * the file spells it: `chr6` on `GRCh38` is `GRCh38#0#chr6`. Undefined when
   * the file holds neither.
   */
  async sequenceFor(
    sample: string,
    refName: string,
    opts: { signal?: AbortSignal } = {},
  ) {
    if (!this.lookupRead) {
      const read = this.walks.getReferenceSequenceNames().then(refNameLookup)
      read.catch(() => {
        if (this.lookupRead === read) {
          this.lookupRead = undefined
        }
      })
      this.lookupRead = read
    }
    const lookup = await raceAbort(this.lookupRead, opts.signal)
    return resolveRefName(lookup, sample, refName)
  }

  async queryStart(start: number, opts: { signal?: AbortSignal } = {}) {
    const { chunk, lookback } = await this.header(opts)
    return chunkQueryStart(start, chunk, lookback)
  }

  async cut(
    refName: string,
    region: { start: number; end: number },
    opts: WalkCutOptions = {},
  ) {
    const {
      signal,
      keep,
      context = WALK_CONTEXT,
      stepBudget = Infinity,
      byteBudget = Infinity,
      links = true,
    } = opts
    const { chunk, lookback } = await this.header({ signal })
    const from = chunkQueryStart(region.start, chunk, lookback)
    const bytesForRegions = this.walks.bytesForRegions?.bind(this.walks)
    const budgeted =
      Number.isFinite(byteBudget) && bytesForRegions !== undefined
    const bytesFor = async (spans: { start: number; end: number }[]) =>
      budgeted
        ? bytesForRegions(
            spans.map(span => ({ refName, ...span })),
            { signal },
          )
        : 0
    const tooLarge = budgeted
      ? await byteBudgetError(
          async end => {
            checkRange(refName, from, end)
            return bytesFor([{ start: from, end }])
          },
          byteBudget,
          region,
          chunk,
          opts.budgetNames,
        )
      : undefined
    if (tooLarge) {
      throw tooLarge
    }
    const reads = new AbortController()
    const stop = () => {
      reads.abort(signal?.reason)
    }
    signal?.throwIfAborted()
    signal?.addEventListener('abort', stop, { once: true })
    const rows: WalkRow[] = []
    let cut: ReturnType<typeof walkCut>
    const graph = new WalkGraph({
      refName,
      start: region.start,
      end: region.end,
      context,
      sequences: opts.sequences,
    })
    const decoded = new Set<number>()
    const ends: DecodedChunk['ends'] = []
    // each chunk read is decoded once, into its rows, nodes and the links
    // between pieces
    const readChunks = async (start: number, end: number) => {
      checkRange(refName, start, end)
      const byChunk = new Map<number, string[][]>()
      await this.walks.getLines(refName, start, end, {
        signal: reads.signal,
        lineCallback: line => {
          const cols = line.split('\t')
          const cs = +cols[1]!
          const held = byChunk.get(cs)
          if (held) {
            held.push(cols)
          } else {
            byChunk.set(cs, [cols])
          }
        },
      })
      for (const [cs, chunkRows] of byChunk) {
        if (decoded.has(cs)) {
          continue
        }
        decoded.add(cs)
        const chunkCut = decodeChunk(chunkRows, keep, opts.sequences)
        rows.push(...chunkCut.rows)
        for (const [id, node] of chunkCut.nodes) {
          graph.placeNode(id, node)
        }
        graph.addSequences(chunkCut.sequences)
        graph.addChunk(chunkCut.chunk)
        ends.push(...chunkCut.ends)
      }
      if (links) {
        for (const [a, b] of pieceLinks(ends)) {
          graph.addLinkEnds(a, b)
        }
      }
    }
    try {
      await readChunks(from, region.end)
      const tooManyRead = stepBudgetError(
        rows,
        stepBudget,
        READ_STEPS_PER_KEPT,
        region,
        chunk,
        keep !== undefined,
        opts.budgetNames,
      )
      if (tooManyRead) {
        throw tooManyRead
      }
      // the chunks read, first to last, which a follow widens to the
      // farthest chunk a cut-off walk goes on in, so the reference between
      // is read too
      let first = from
      let last = Math.floor((region.end - 1) / chunk) * chunk
      const reach = FOLLOW_CHUNKS * chunk
      const [lo, hi] = [first - reach, last + reach]
      cut = walkCut(joinPieces(rows), graph.nodes, refName)
      let spent = await bytesFor([{ start: from, end: region.end }])
      for (let round = 0; round < FOLLOW_ROUNDS; round++) {
        const starts = cut.open
          .filter(c => c.refName === refName && c.start >= lo && c.start <= hi)
          .map(c => c.start)
        const spans = [
          { start: Math.min(first, ...starts), end: first },
          { start: last + chunk, end: Math.max(last, ...starts) + 1 },
        ].filter(span => span.end > span.start)
        if (spans.length === 0) {
          break
        }
        spent += await bytesFor(spans)
        if (spent > byteBudget) {
          break
        }
        first = Math.min(first, ...starts)
        last = Math.max(last, ...starts)
        await Promise.all(spans.map(({ start, end }) => readChunks(start, end)))
        cut = walkCut(joinPieces(rows), graph.nodes, refName)
      }
    } catch (error) {
      reads.abort()
      throw error
    } finally {
      signal?.removeEventListener('abort', stop)
    }
    const tooMany = keptStepsError(
      cut.fragments,
      stepBudget,
      region,
      keep !== undefined,
    )
    if (tooMany) {
      throw tooMany
    }
    return { graph, ...cut }
  }

  /**
   * The graph under a region as tables: the reference sample's sequence over
   * `region`, with `context` bp either side, and the walks of the haplotypes
   * asked for. `region.refName` may be a contig (`chr6`) or the file's whole
   * name; the result names the sequence as the file does.
   */
  async cutRegion(
    region: { refName: string; start: number; end: number },
    opts: RegionCutOptions = {},
  ) {
    const { haplotypes, referenceSample, nodeLimit = Infinity, ...rest } = opts
    const { signal } = opts
    const header = await this.header({ signal })
    const sample = referenceSample ?? header.references[0] ?? ''
    const refName = await this.sequenceFor(sample, region.refName, { signal })
    if (refName === undefined) {
      const names = await this.referenceSequences({ signal })
      throw new Error(
        `The walk file has no sequence ${region.refName} on ${sample}; it holds ${names.slice(0, 5).join(', ')}${names.length > 5 ? ', …' : ''}`,
      )
    }
    const { graph, kept, fragments } = await this.cut(refName, region, {
      ...rest,
      keep: walksKept(haplotypes, refName),
    })
    const tooMany = keptNodesError(
      kept.size,
      nodeLimit,
      region,
      !!haplotypes?.length,
    )
    if (tooMany) {
      throw tooMany
    }
    return { tables: graph.tables(kept, fragments), refName, sample }
  }
}
