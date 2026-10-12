import { decodeChunk, pieceLinks } from './walkFormat2.ts'
import {
  WalkGraph,
  byteBudgetError,
  chunkQueryStart,
  joinPieces,
  keptStepsError,
  lookbackChunks,
  parseWalkRow,
  stepBudgetError,
  walkCut,
  walkHeader,
  walkRowName,
} from './walkRows.ts'

import type { DecodedChunk } from './walkFormat2.ts'
import type { WalkHeader, WalkRow } from './walkRows.ts'

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

// Format 1 is three files, `<prefix>.<sample>.walks/nodes/links.bed.gz`;
// format 2 is the walk file alone
export interface WalkFileSet {
  walks: TabixSource
  nodes?: TabixSource
  links?: TabixSource
}

// @gmod/tabix loops forever on a NaN end and reads nothing for a NaN start
// (GMOD/tabix-js#157), so every query refuses a range that is not finite first
function checkRange(refName: string, start: number, end: number) {
  if (!Number.isFinite(start) || !Number.isFinite(end)) {
    throw new Error(
      `${refName}:${start}-${end} is not a finite range to query the graph's index for`,
    )
  }
}

function getLines(
  file: TabixSource,
  refName: string,
  start: number,
  end: number,
  opts: { signal?: AbortSignal; lineCallback: (line: string) => void },
) {
  checkRange(refName, start, end)
  return file.getLines(refName, start, end, opts)
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
  // false leaves the link file unread, for a caller that draws no edges
  links?: boolean
  // keep the bases of format 1 nodes written with --sequences
  sequences?: boolean
}

export type WalkReaderHeader = WalkHeader & { chunk: number; lookback: number }

/**
 * A walk-indexed graph's three files, read a window at a time: one read of
 * each over whole chunks, from the lookback chunks before the window to its
 * end. The indexes first estimate the bytes those reads fetch, and a window
 * past the byte budget reads no row. Only the walks `keep` accepts are
 * decoded. The step budget bounds the steps the cut keeps, and rows holding
 * READ_STEPS_PER_KEPT times as many stop the reads, which go out together.
 */
export class WalkReader {
  private headerRead: Promise<WalkReaderHeader> | undefined

  // `chunkWithoutHeader` is the chunk size of a walk file from before
  // gfa-to-tabix 0.4.0, whose header does not state one
  constructor(
    private files: WalkFileSet,
    private chunkWithoutHeader = 65536,
  ) {}

  // The header, read once for every caller: one caller giving up stops only
  // its own wait, and a read that fails is tried again by the next caller
  header(opts: { signal?: AbortSignal } = {}) {
    if (!this.headerRead) {
      const read = this.files.walks.getHeaderLines().then(lines => {
        const header = walkHeader(lines)
        const chunk = header.chunk ?? this.chunkWithoutHeader
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
    return this.files.walks.getReferenceSequenceNames(opts)
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
    } = opts
    const { chunk, lookback, format } = await this.header({ signal })
    // format 2 holds the whole cut in its walk file
    const format2 = format === 2
    const nodeFile = this.files.nodes
    const linkFile = opts.links === false ? undefined : this.files.links
    if (!format2 && !nodeFile) {
      throw new Error(
        'A format 1 walk file needs its node file: open it with `nodes`',
      )
    }
    const links = format2 ? opts.links !== false : linkFile !== undefined
    const from = chunkQueryStart(region.start, chunk, lookback)
    const files = format2
      ? [this.files.walks]
      : [this.files.walks, nodeFile!, ...(linkFile ? [linkFile] : [])]
    const budgeted =
      Number.isFinite(byteBudget) &&
      files.every(file => file.bytesForRegions !== undefined)
    const bytesFor = async (spans: { start: number; end: number }[]) => {
      if (!budgeted) {
        return 0
      }
      const sizes = await Promise.all(
        files.map(file =>
          file.bytesForRegions!(
            spans.map(span => ({ refName, ...span })),
            { signal },
          ),
        ),
      )
      return sizes.reduce((sum, bytes) => sum + bytes, 0)
    }
    const tooLarge = budgeted
      ? await byteBudgetError(
          async end => {
            checkRange(refName, from, end)
            return bytesFor([{ start: from, end }])
          },
          byteBudget,
          region,
          chunk,
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
    const read = (
      file: TabixSource,
      lineCallback: (line: string) => void,
      start = from,
      end = region.end,
    ) =>
      getLines(file, refName, start, end, {
        signal: reads.signal,
        lineCallback,
      })
    const addRow = (line: string) => {
      if (keep === undefined || keep(walkRowName(line))) {
        rows.push(parseWalkRow(line))
      }
    }
    const addNode = (line: string) => {
      graph.addNode(line)
    }
    const addLink = (line: string) => {
      graph.addLink(line)
    }
    const decoded = new Set<number>()
    const ends: DecodedChunk['ends'] = []
    // Format 2: each chunk read is decoded once, into the rows, nodes and
    // links format 1 would have filed, and the links between pieces read
    const readChunks = async (start = from, end = region.end) => {
      const byChunk = new Map<number, string[][]>()
      await read(
        this.files.walks,
        line => {
          const cols = line.split('\t')
          const cs = +cols[1]!
          const held = byChunk.get(cs)
          if (held) {
            held.push(cols)
          } else {
            byChunk.set(cs, [cols])
          }
        },
        start,
        end,
      )
      for (const [cs, chunkRows] of byChunk) {
        if (decoded.has(cs)) {
          continue
        }
        decoded.add(cs)
        const chunkCut = decodeChunk(chunkRows, keep)
        rows.push(...chunkCut.rows)
        for (const [id, node] of chunkCut.nodes) {
          graph.placeNode(id, node)
        }
        graph.addFormat2Chunk(chunkCut.chunk)
        ends.push(...chunkCut.ends)
      }
      if (links) {
        for (const [a, b] of pieceLinks(ends)) {
          graph.addLinkEnds(a, b)
        }
      }
    }
    const checkSteps = () => {
      const error = stepBudgetError(
        rows,
        stepBudget,
        READ_STEPS_PER_KEPT,
        region,
        chunk,
        keep !== undefined,
      )
      if (error) {
        throw error
      }
    }
    try {
      await (format2
        ? readChunks().then(checkSteps)
        : Promise.all([
            read(this.files.walks, addRow).then(checkSteps),
            read(nodeFile!, addNode),
            ...(linkFile ? [read(linkFile, addLink)] : []),
          ]))
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
        await Promise.all(
          spans.flatMap(({ start, end }) =>
            format2
              ? [readChunks(start, end)]
              : [
                  read(this.files.walks, addRow, start, end),
                  read(nodeFile!, addNode, start, end),
                  ...(linkFile ? [read(linkFile, addLink, start, end)] : []),
                ],
          ),
        )
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
}
