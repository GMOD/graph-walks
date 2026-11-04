import { refNameLookup, resolveRefName } from './refNames.ts'
import { chunkRef, decodeChunk, pieceLinks } from './walkChunk.ts'
import {
  ATTACH_BP,
  departureRowKey,
  parseDepartureRow,
  squeezer,
} from './walkDeparture.ts'
import {
  WalkGraph,
  byteBudgetError,
  chunkQueryStart,
  joinPieces,
  keptNodesError,
  keptStepsError,
  lookbackChunks,
  megabytes,
  stepBudgetError,
  surelyKept,
  walkCut,
  walkHeader,
  walkQueryStart,
  walksKept,
} from './walkRows.ts'

import type { ChunkRef, DecodedChunk } from './walkChunk.ts'
import type { WalkDeparture } from './walkDeparture.ts'
import type { BudgetNames, WalkHeader, WalkRow } from './walkRows.ts'

// A window reads whole chunks and keeps the steps on it, so a narrow window
// reads several times what it keeps: 464 haplotypes across the 1.5 kb ABCA7
// repeat read 5.0 M steps and keep 1.16 M, and the cut takes 0.4 s from warm
// files. The step budget bounds what the cut keeps, and rows are read to this
// multiple of it.
const READ_STEPS_PER_KEPT = 8

// bp of reference a cut keeps past each end of the window
const WALK_CONTEXT = 1000

// detail records an expanding cut reads at once
const DETAIL_READS = 16

// chunks a cut fetches for the own steps of departure records whose leave
// piece its read lacks
const ATTACH_FETCHES = 8

// The budgets the JBrowse graph plugin and BandageJS cut a browser's window by
export const WALK_BUDGETS = {
  stepBudget: 4_000_000,
  byteBudget: 8_000_000,
  departureByteBudget: 2_000_000,
  attachBp: ATTACH_BP,
}

/**
 * A Tabix-indexed file as the reader uses one: @gmod/tabix's
 * TabixIndexedFile is one. `bytesForRegions` estimates the compressed bytes a
 * read fetches, which the byte budgets need.
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

// A walk file in a format other than 6, which the reader refuses
export class WalkFormatError extends Error {
  override name = 'WalkFormatError'
}

/**
 * A window whose departure read would fetch more than its budget. Zooming in
 * does not lower the read, since a record is filed over its whole stretch,
 * so this is not a RegionTooLargeError and shows no zoom-in notice.
 */
export class DepartureReadError extends Error {
  override name = 'DepartureReadError'
  readonly bytes: number
  readonly limit: number

  constructor(bytes: number, limit: number) {
    super(
      `Too many departures from the reference over this window to read (${megabytes(bytes)} against a budget of ${megabytes(limit)}). A collapsed repeat usually causes this: haplotypes hop between its copies, and each hop is a departure`,
    )
    this.bytes = bytes
    this.limit = limit
  }
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

// A detail record's sequence is its chunk's, `~` and the chunk's start, and
// a sequence's departures lie under `~departure`
function withoutSiblings(names: string[]) {
  return names.filter(name => !name.includes('~'))
}

// Reads `items`, at most `limit` at once and none once `signal` aborts, and
// hands each to `use` in order
export async function readInOrder<T, R>(
  items: T[],
  limit: number,
  signal: AbortSignal,
  read: (item: T) => Promise<R>,
  use: (item: T, result: R) => void,
) {
  const reads: Promise<R>[] = []
  const next = () => {
    if (reads.length < items.length && !signal.aborted) {
      const started = read(items[reads.length]!)
      reads.push(started)
      started.then(next, () => {})
    }
  }
  for (let k = 0; k < limit; k++) {
    next()
  }
  for (let i = 0; i < items.length; i++) {
    signal.throwIfAborted()
    use(items[i]!, await reads[i]!)
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
  // the most compressed bytes a cut fetches from the windows' chunks;
  // Infinity, the default, has no limit. The detail records `expand` reads
  // count only against stepBudget
  byteBudget?: number
  // the most compressed bytes the departure read fetches; Infinity, the
  // default, has no limit
  departureByteBudget?: number
  // the most bp of own sequence a departure gets its steps attached for,
  // ATTACH_BP by default
  attachBp?: number
  // the most chunks the cut fetches for the own steps of records whose leave
  // piece its read lacks
  attachFetches?: number
  context?: number
  // false lists no links, for a caller that draws no edges
  links?: boolean
  // keep the bases of the nodes, from a file written with --sequences
  sequences?: boolean
  // how a zoom-in notice names the budgets
  budgetNames?: BudgetNames
  // the most nodes a cut keeps; Infinity, the default, has no limit
  nodeLimit?: number
  // keep the reference over the window and its context alone, and split each
  // walk where it steps onto the reference off them
  contained?: boolean
  // the abridged walks to read in full from their detail records, by
  // `sample#haplotype#contig`; none by default
  expand?: (name: string) => boolean
}

export interface RegionCutOptions extends Omit<WalkCutOptions, 'keep'> {
  // the `sample#haplotype` walks to keep beside the reference; every walk when
  // undefined or empty
  haplotypes?: string[]
  // the sample whose sequence the region names, the file's first reference
  // by default
  referenceSample?: string
}

export type WalkReaderHeader = WalkHeader & { chunk: number }

/**
 * A `gfa-to-tabix --walks` file, `<prefix>.<sample>.walks.bed.gz`, read a
 * window at a time: the walk read over whole chunks, from the lookback
 * chunks before the window to its end, each chunk decoded once, and the
 * departure read over the window and its context on the sequence's
 * `~departure` sibling. The index first estimates the bytes each fetches,
 * and a window past a budget reads no row. Only the walks `keep` accepts are
 * decoded. The step budget bounds the steps the cut keeps, and rows holding
 * READ_STEPS_PER_KEPT times as many stop the read.
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
        if (header.format !== 6) {
          const format = header.format ?? 1
          throw new WalkFormatError(
            format > 6
              ? `This walk file is format ${format}; this version of @jbrowse/graph-walks reads format 6, so update the viewer, or @jbrowse/graph-walks, to one that reads format ${format}`
              : `This walk file is format ${format}; @jbrowse/graph-walks reads format 6, so rebuild it with gfa-to-tabix 3.0 or later`,
          )
        }
        return { ...header, chunk: header.chunk ?? 65536 }
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

  // The reference sequences the walk file holds windows on, without the
  // `seq~cs` sequences of detail records and the `seq~departure` siblings
  async referenceSequences(opts: { signal?: AbortSignal } = {}) {
    return withoutSiblings(await this.walks.getReferenceSequenceNames(opts))
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
      const read = this.walks
        .getReferenceSequenceNames()
        .then(names => refNameLookup(withoutSiblings(names)))
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

  // The first base a cut over `start`, and `end` on `refName` when given,
  // queries
  async queryStart(
    start: number,
    opts: { signal?: AbortSignal; refName?: string; end?: number } = {},
  ) {
    const { refName = '', end = start + 1 } = opts
    return walkQueryStart(await this.header(opts), refName, start, end)
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
      departureByteBudget = Infinity,
      attachBp = ATTACH_BP,
      attachFetches = ATTACH_FETCHES,
      links = true,
      contained = false,
      nodeLimit = Infinity,
      expand,
    } = opts
    const header = await this.header({ signal })
    const { chunk } = header
    const from = walkQueryStart(header, refName, region.start, region.end)
    const sibling = `${refName}~departure`
    // the drawn stretch reaches at most a node past the window's context
    const reach = header.maxNode ?? chunk
    const extent = {
      start: Math.max(0, region.start - context - reach),
      end: region.end + context + reach,
    }
    const bytesForRegions = this.walks.bytesForRegions?.bind(this.walks)
    const bytesFor = async (
      name: string,
      spans: { start: number; end: number }[],
    ) =>
      bytesForRegions === undefined
        ? 0
        : bytesForRegions(
            spans.map(span => ({ refName: name, ...span })),
            { signal },
          )
    if (Number.isFinite(byteBudget) && bytesForRegions !== undefined) {
      const tooLarge = await byteBudgetError(
        async end => {
          checkRange(refName, from, end)
          return bytesFor(refName, [{ start: from, end }])
        },
        byteBudget,
        region,
        chunk,
        opts.budgetNames,
      )
      if (tooLarge) {
        throw tooLarge
      }
    }
    if (Number.isFinite(departureByteBudget) && bytesForRegions !== undefined) {
      checkRange(sibling, extent.start, extent.end)
      const bytes = await bytesFor(sibling, [extent])
      if (bytes > departureByteBudget) {
        throw new DepartureReadError(bytes, departureByteBudget)
      }
    }
    const reads = new AbortController()
    const stop = () => {
      reads.abort(signal?.reason)
    }
    signal?.throwIfAborted()
    signal?.addEventListener('abort', stop, { once: true })
    const rows: WalkRow[] = []
    const graph = new WalkGraph({
      refName,
      start: region.start,
      end: region.end,
      context,
      sequences: opts.sequences,
    })
    const decoded = new Set<number>()
    const ends: DecodedChunk['ends'] = []
    const surely = new Set<number>()
    const bounds = {
      contained,
      start: region.start - context,
      end: region.end + context,
    }
    const filtered = keep !== undefined
    // the expanded walks' rows, with the R of their chunk
    let details: { ref: ChunkRef; row: WalkRow }[] = []
    // each chunk read is decoded once, into its rows, nodes and the links
    // between pieces, and a cut over its node limit stops there
    const decode = (
      cs: number,
      chunkRows: string[][],
      only: ((name: string) => boolean) | undefined,
    ) => {
      if (decoded.has(cs)) {
        return
      }
      decoded.add(cs)
      const chunkCut = decodeChunk(chunkRows, only ?? keep, opts.sequences)
      let ref: ChunkRef | undefined
      for (const row of chunkCut.rows) {
        if (row.abridged && expand?.(row.name)) {
          ref ??= chunkRef(chunkRows)
          details.push({ ref, row })
        } else {
          rows.push(row)
        }
      }
      place(chunkCut, cs)
    }
    const readDetail = async ({ row }: { row: WalkRow }) => {
      const { refName: sequence, start } = row.abridged!.detail
      checkRange(sequence, start, start + 1)
      const held: string[][] = []
      await this.walks.getLines(sequence, start, start + 1, {
        signal: reads.signal,
        lineCallback: line => {
          const cols = line.split('\t')
          cols[0] = refName
          cols[1] = String(row.chunkStart)
          held.push(cols)
        },
      })
      return held
    }
    const decodeDetail = (
      { ref, row }: { ref: ChunkRef; row: WalkRow },
      held: string[][],
    ) => {
      const chunkCut = decodeChunk(
        held,
        name => name === row.name,
        opts.sequences,
        ref,
      )
      for (const decodedRow of chunkCut.rows) {
        if (decodedRow.name === row.name) {
          rows.push(decodedRow)
        }
      }
      place(chunkCut, row.chunkStart)
    }
    const place = (chunkCut: DecodedChunk, cs: number) => {
      for (const [id, node] of chunkCut.nodes) {
        graph.placeNode(id, node)
      }
      for (const { id, ...node } of chunkCut.revisits) {
        graph.placeNode(id, node)
      }
      graph.addSequences(chunkCut.sequences)
      graph.addChunk(chunkCut.chunk)
      for (const end of chunkCut.ends) {
        ends.push(end)
      }
      if (Number.isFinite(nodeLimit)) {
        surelyKept(chunkCut, graph.nodes, bounds, surely)
        if (surely.size > nodeLimit) {
          surely.clear()
          surelyKept({ nodes: graph.nodes, rows }, graph.nodes, bounds, surely)
        }
        const over = keptNodesError(
          surely.size,
          nodeLimit,
          region,
          filtered,
          cs + chunk,
        )
        if (over) {
          throw over
        }
      }
    }
    // rows arrive in chunk order, and each chunk is decoded as the next
    // starts; `only` narrows a fetch to the walks named
    const readChunks = async (
      start: number,
      end: number,
      only?: (name: string) => boolean,
    ) => {
      checkRange(refName, start, end)
      let at: number | undefined
      let held: string[][] = []
      await this.walks.getLines(refName, start, end, {
        signal: reads.signal,
        lineCallback: line => {
          const cols = line.split('\t')
          const cs = +cols[1]!
          if (cs !== at) {
            if (at !== undefined) {
              decode(at, held, only)
            }
            at = cs
            held = []
          }
          held.push(cols)
        },
      })
      if (at !== undefined) {
        decode(at, held, only)
      }
      const pending = details
      details = []
      await readInOrder(
        pending,
        DETAIL_READS,
        reads.signal,
        readDetail,
        decodeDetail,
      )
      if (links) {
        for (const [a, b] of pieceLinks(ends)) {
          graph.addLinkEnds(a, b)
        }
      }
    }
    // the records over the extent, each once, of the walks kept
    const readDepartures = async () => {
      checkRange(sibling, extent.start, extent.end)
      const seen = new Set<string>()
      const records: WalkDeparture[] = []
      await this.walks.getLines(sibling, extent.start, extent.end, {
        signal: reads.signal,
        lineCallback: line => {
          const cols = line.split('\t')
          if (cols[3] !== 'departure') {
            return
          }
          const key = departureRowKey(cols)
          if (seen.has(key)) {
            return
          }
          seen.add(key)
          const record = parseDepartureRow(cols)
          if (keep === undefined || keep(record.name)) {
            records.push(record)
          }
        },
      })
      return records
    }
    let cut: ReturnType<typeof walkCut>
    try {
      const [records] = await Promise.all([
        readDepartures(),
        readChunks(from, region.end),
      ])
      const tooManyRead = stepBudgetError(
        rows,
        stepBudget,
        READ_STEPS_PER_KEPT,
        region,
        chunk,
        filtered,
        opts.budgetNames,
      )
      if (tooManyRead) {
        throw tooManyRead
      }
      const cutRows = () =>
        walkCut(joinPieces(rows), graph.nodes, refName, records, {
          contained,
          context,
          attachBp,
          squeeze: squeezer(header.holes, refName),
        })
      cut = cutRows()
      // the leave pieces the read lacks, each chunk fetched once for the
      // records wanting it, with the chunk before for a straddling node
      const fetches = new Map<number, Set<string>>()
      for (const { departure, base } of cut.unattached) {
        const cs = Math.floor(base / chunk) * chunk
        if (!decoded.has(cs)) {
          fetches.set(cs, (fetches.get(cs) ?? new Set()).add(departure.name))
        }
      }
      const wanted = [...fetches].slice(0, attachFetches)
      if (wanted.length > 0) {
        await Promise.all(
          wanted.map(([cs, names]) =>
            readChunks(
              chunkQueryStart(cs, chunk, lookbackChunks(header.maxNode, chunk)),
              cs + 1,
              name => names.has(name),
            ),
          ),
        )
        cut = cutRows()
      }
    } catch (error) {
      reads.abort()
      throw error
    } finally {
      signal?.removeEventListener('abort', stop)
    }
    const tooMany =
      keptStepsError(
        cut.fragments,
        cut.departures,
        stepBudget,
        region,
        filtered,
      ) ?? keptNodesError(cut.kept.size, nodeLimit, region, filtered)
    if (tooMany) {
      throw tooMany
    }
    const { kept, fragments, departures } = cut
    return { graph, kept, fragments, departures }
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
    const { haplotypes, referenceSample, ...rest } = opts
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
    const { graph, kept, fragments, departures } = await this.cut(
      refName,
      region,
      { ...rest, keep: walksKept(haplotypes, refName) },
    )
    return {
      tables: graph.tables(kept, fragments, departures),
      refName,
      sample,
    }
  }
}
