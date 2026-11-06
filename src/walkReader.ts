import { RemoteFile } from 'generic-filehandle2'

import { refNameLookup, resolveRefName } from './refNames.ts'
import { checkRange } from './tabixSource.ts'
import {
  chunkRef,
  decodeChunk,
  pieceLinks,
  referenceNodes,
} from './walkChunk.ts'
import {
  ATTACH_BP,
  departureOwnTiles,
  parseDepartureRow,
  squeezer,
} from './walkDeparture.ts'
import { WalkFile, WalkFormatError, once, raceAbort } from './walkFile.ts'
import {
  WalkGraph,
  byteBudgetError,
  joinPieces,
  keptNodesError,
  keptStepsError,
  megabytes,
  stepBudgetError,
  surelyKept,
  walkCut,
  walkHeader,
  walksKept,
} from './walkRows.ts'

import type { ChunkRef, DecodedChunk } from './walkChunk.ts'
import type { WalkDeparture } from './walkDeparture.ts'
import type { WalkFilehandle, WalkUnit } from './walkFile.ts'
import type { BudgetNames, WalkHeader, WalkRow } from './walkRows.ts'

// A window reads whole tiles and keeps the steps on it, so a narrow window
// reads several times what it keeps: 464 haplotypes across the 1.5 kb ABCA7
// repeat read 5.0 M steps and keep 1.16 M, and the cut takes 0.4 s from warm
// files. The step budget bounds what the cut keeps, and rows are read to this
// multiple of it.
const READ_STEPS_PER_KEPT = 8

// bp of reference a cut keeps past each end of the window
const WALK_CONTEXT = 1000

// detail records an expanding cut reads at once
const DETAIL_READS = 16

// tiles a cut fetches for the own steps of departure records whose leave
// piece its read lacks
const ATTACH_FETCHES = 8

// The budgets a cut keeps to unless its caller passes others: what a browser
// window can hold. Infinity lifts one
export const WALK_BUDGETS = {
  stepBudget: 4_000_000,
  byteBudget: 8_000_000,
  departureByteBudget: 2_000_000,
  attachBp: ATTACH_BP,
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

export interface WalkCutOptions {
  signal?: AbortSignal
  // which walks to decode, by `sample#haplotype#contig`; undefined is all
  keep?: (name: string) => boolean
  // the most steps a cut keeps, WALK_BUDGETS' by default; Infinity has no
  // limit
  stepBudget?: number
  // the most compressed bytes a cut fetches from the tiles under its window
  // and context, WALK_BUDGETS' by default; Infinity has no limit. The detail
  // records `expand` reads count only against stepBudget, and the tiles
  // fetched for records' own steps against no budget
  byteBudget?: number
  // the most compressed bytes the departure read fetches, WALK_BUDGETS' by
  // default; Infinity has no limit
  departureByteBudget?: number
  // the most bp of own sequence a departure gets its steps attached for,
  // ATTACH_BP by default
  attachBp?: number
  // the most tiles the cut fetches for the own steps of records whose leave
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

export type WalkReaderHeader = WalkHeader & { tile: number }

const overlaps = (unit: WalkUnit, span: { start: number; end: number }) =>
  unit.length > 0 && unit.start < span.end && unit.end > span.start

const bytesOf = (units: WalkUnit[]) =>
  units.reduce((sum, unit) => sum + unit.length, 0)

// The stretch a cut draws: the window and its context
const drawnSpan = (
  region: { start: number; end: number },
  context: number,
) => ({
  start: Math.max(0, region.start - context),
  end: region.end + context,
})

// The tiles an attach fetch reads for a record: those its `ownTiles` names
export function attachTiles(units: WalkUnit[], ownTiles: number[]) {
  return units.filter(
    u => u.kind === 'tile' && u.length > 0 && ownTiles.includes(Number(u.name)),
  )
}

// The directory beside a data file's URL, before any query or fragment
const directoryOf = (url: string) => url.replace(/(?=$|[?#])/, '.idx')

// The records of departure lines over `span`, each line once, of the walks
// kept, with the tiles holding each one's own steps
function departureRecords(
  lines: string[],
  span: { start: number; end: number },
  keep: ((name: string) => boolean) | undefined,
) {
  const records: WalkDeparture[] = []
  const ownTiles = new Map<WalkDeparture, number[]>()
  for (const line of new Set(lines)) {
    const cols = line.split('\t')
    if (
      cols[3] === 'departure' &&
      +cols[1]! < span.end &&
      +cols[2]! > span.start
    ) {
      const record = parseDepartureRow(cols)
      if (keep === undefined || keep(record.name)) {
        records.push(record)
        ownTiles.set(record, departureOwnTiles(cols))
      }
    }
  }
  return { records, ownTiles }
}

/**
 * A `gfa-to-tabix --walks` file of walk format 7, `<prefix>.<sample>.walks.gz`
 * and its directory `<prefix>.<sample>.walks.gz.idx`, read a window at a
 * time. At open the directory, the header and the sequence's unit table; for
 * a window, the tiles and departure blocks overlapping the window and its
 * context, in one read wherever their bytes touch, each tile decoded once.
 * The table gives each read's bytes before it is made, and a window past a
 * budget reads no row. Only the walks `keep` accepts are decoded. The step
 * budget bounds the steps the cut keeps, and rows holding
 * READ_STEPS_PER_KEPT times as many stop the read.
 */
export class WalkReader {
  private readonly file: WalkFile

  constructor(file: WalkFile) {
    this.file = file
  }

  /**
   * The walk file at `url` and its directory at `url` plus `.idx`, before any
   * query string, each opened by `open`, a RemoteFile by default
   */
  static open(
    url: string,
    open: (url: string) => WalkFilehandle = location =>
      new RemoteFile(location),
  ) {
    return WalkReader.fromFilehandles({
      data: open(url),
      directory: open(directoryOf(url)),
    })
  }

  // The walk file `data` with its directory, the `.idx` file
  static fromFilehandles({
    data,
    directory,
  }: {
    data: WalkFilehandle
    directory: WalkFilehandle
  }) {
    return new WalkReader(new WalkFile(data, directory))
  }

  private readonly readHeader = once(async () => {
    const header = walkHeader(await this.file.header())
    const { tile } = header
    if (tile === undefined) {
      throw new WalkFormatError("The walk file's #walks line states no tile:i:")
    }
    return { ...header, tile }
  })

  private readonly readLookup = once(async () =>
    refNameLookup(await this.file.sequences()),
  )

  // The header, read once for every caller: one caller giving up stops only
  // its own wait, and a read that fails is tried again by the next caller
  header(opts: { signal?: AbortSignal } = {}) {
    return raceAbort(this.readHeader(), opts.signal)
  }

  // The reference sequences the walk file holds windows on
  referenceSequences(opts: { signal?: AbortSignal } = {}) {
    return raceAbort(this.file.sequences(), opts.signal)
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
    const lookup = await raceAbort(this.readLookup(), opts.signal)
    return resolveRefName(lookup, sample, refName)
  }

  // The units the file holds for a reference sequence, in file order
  units(refName: string, opts: { signal?: AbortSignal } = {}) {
    return raceAbort(this.file.units(refName), opts.signal)
  }

  // Where the reference's last node on `refName` ends, from the extents of
  // its tiles; 0 for a sequence the file lacks
  async sequenceLength(refName: string, opts: { signal?: AbortSignal } = {}) {
    const units = await this.units(refName, opts)
    return units.reduce(
      (end, u) => (u.kind === 'tile' ? Math.max(end, u.end) : end),
      0,
    )
  }

  // The compressed bytes a cut of `region` reads before any attach fetch or
  // expansion: its tiles and departure blocks
  async bytesFor(
    refName: string,
    region: { start: number; end: number },
    opts: { signal?: AbortSignal; context?: number } = {},
  ) {
    const drawn = drawnSpan(region, opts.context ?? WALK_CONTEXT)
    const units = await this.units(refName, opts)
    return bytesOf(
      units.filter(
        u =>
          (u.kind === 'tile' || u.kind === 'departures') && overlaps(u, drawn),
      ),
    )
  }

  // Every pass of the reference over `region`, by start, from the ref rows of
  // the tiles there: a node the reference passes twice there comes twice
  async referenceNodes(
    refName: string,
    region: { start: number; end: number },
    opts: { signal?: AbortSignal } = {},
  ) {
    checkRange(refName, region.start, region.end)
    const tiles = (await this.units(refName, opts)).filter(
      u => u.kind === 'tile' && overlaps(u, region),
    )
    const read = await this.file.read(tiles, opts.signal)
    return referenceNodes(read.flat().map(line => line.split('\t')))
      .filter(
        node =>
          node.refName === refName &&
          node.end > region.start &&
          node.start < region.end,
      )
      .sort((a, b) => a.start - b.start || a.end - b.end)
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
      stepBudget = WALK_BUDGETS.stepBudget,
      byteBudget = WALK_BUDGETS.byteBudget,
      departureByteBudget = WALK_BUDGETS.departureByteBudget,
      attachBp = ATTACH_BP,
      attachFetches = ATTACH_FETCHES,
      links = true,
      contained = false,
      nodeLimit = Infinity,
      expand,
    } = opts
    checkRange(refName, region.start, region.end)
    const [header, units] = await Promise.all([
      this.header({ signal }),
      this.units(refName, { signal }),
    ])
    const { tile } = header
    const drawn = drawnSpan(region, context)
    const tiles = units.filter(u => u.kind === 'tile' && overlaps(u, drawn))
    const blocks = units.filter(
      u => u.kind === 'departures' && overlaps(u, drawn),
    )
    if (Number.isFinite(byteBudget)) {
      const tooLarge = byteBudgetError(
        end => bytesOf(tiles.filter(u => u.start < end)),
        byteBudget,
        region,
        tile,
        opts.budgetNames,
      )
      if (tooLarge) {
        throw tooLarge
      }
    }
    const departureBytes = bytesOf(blocks)
    if (departureBytes > departureByteBudget) {
      throw new DepartureReadError(departureBytes, departureByteBudget)
    }
    const reads = new AbortController()
    const stop = () => {
      reads.abort(signal?.reason)
    }
    signal?.throwIfAborted()
    signal?.addEventListener('abort', stop, { once: true })
    const rows: WalkRow[] = []
    // the rows of tiles fetched to attach, which give records their own steps
    // and nothing else
    const attachRows: WalkRow[] = []
    const graph = new WalkGraph({
      refName,
      start: region.start,
      end: region.end,
      context,
      sequences: opts.sequences,
    })
    // the name filters each tile was decoded under, by its offset
    const decodedBy = new Map<number, ((name: string) => boolean)[]>()
    const decodedFor = (unit: WalkUnit, name: string) =>
      decodedBy.get(unit.offset)?.some(took => took(name)) === true
    const ends: DecodedChunk['ends'] = []
    const surely = new Set<number>()
    const bounds = {
      contained,
      start: region.start - context,
      end: region.end + context,
    }
    const filtered = keep !== undefined
    const detailUnits = new Map(
      units.filter(u => u.kind === 'detail').map(u => [u.name, u]),
    )
    // the detail records of the expanded walks, each once, with the R of
    // their tile and the rows they join; a path's abridged pieces in a tile
    // share one
    let details: { ref: ChunkRef; row: WalkRow; into: WalkRow[] }[] = []
    const detailsRead = new Set<string>()
    const place = (tileCut: DecodedChunk, readTo: number) => {
      for (const [id, node] of tileCut.nodes) {
        graph.placeNode(id, node)
      }
      for (const { id, ...node } of tileCut.revisits) {
        graph.placeNode(id, node)
      }
      graph.addSequences(tileCut.sequences)
      graph.addChunk(tileCut.chunk)
      for (const end of tileCut.ends) {
        ends.push(end)
      }
      if (Number.isFinite(nodeLimit)) {
        surelyKept(tileCut, graph.nodes, bounds, surely)
        if (surely.size > nodeLimit) {
          surely.clear()
          surelyKept({ nodes: graph.nodes, rows }, graph.nodes, bounds, surely)
        }
        const over = keptNodesError(
          surely.size,
          nodeLimit,
          region,
          filtered,
          readTo,
        )
        if (over) {
          throw over
        }
      }
    }
    // a tile is decoded into `into`, with its nodes and the links between
    // pieces, for the names `only` takes that no earlier decode of it took,
    // and a cut over its node limit stops there
    const decode = (
      unit: WalkUnit,
      lines: string[],
      into: WalkRow[],
      only?: (name: string) => boolean,
    ) => {
      const wanted = only ?? keep ?? (() => true)
      const earlier = decodedBy.get(unit.offset) ?? []
      decodedBy.set(unit.offset, [...earlier, wanted])
      const tileRows = lines.map(line => line.split('\t'))
      const tileCut = decodeChunk(
        tileRows,
        Number(unit.name),
        name => wanted(name) && !earlier.some(took => took(name)),
        opts.sequences,
      )
      let ref: ChunkRef | undefined
      for (const row of tileCut.rows) {
        if (row.abridged && expand?.(row.name)) {
          const { name } = row.abridged.detail
          if (!detailsRead.has(name)) {
            detailsRead.add(name)
            ref ??= chunkRef(tileRows)
            details.push({ ref, row, into })
          }
        } else {
          into.push(row)
        }
      }
      place(tileCut, unit.end)
    }
    const readDetail = async ({ row }: { row: WalkRow }) => {
      const { name } = row.abridged!.detail
      const unit = detailUnits.get(name)
      if (!unit) {
        throw new WalkFormatError(
          `${refName}'s unit table lists no detail unit ${name}, which an abridged row of ${row.name} names`,
        )
      }
      const [lines = []] = await this.file.read([unit], reads.signal)
      return lines.map(line => line.split('\t'))
    }
    const decodeDetail = (
      { ref, row, into }: { ref: ChunkRef; row: WalkRow; into: WalkRow[] },
      held: string[][],
    ) => {
      const tileCut = decodeChunk(
        held,
        row.chunkStart,
        name => name === row.name,
        opts.sequences,
        ref,
      )
      for (const decodedRow of tileCut.rows) {
        if (decodedRow.name === row.name) {
          into.push(decodedRow)
        }
      }
      place(tileCut, row.chunkStart + tile)
    }
    // the detail records the tiles decoded asked for, then the links
    // between pieces
    const finish = async () => {
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
    let cut: ReturnType<typeof walkCut>
    try {
      const read = await this.file.read([...tiles, ...blocks], reads.signal)
      tiles.forEach((unit, k) => {
        decode(unit, read[k]!, rows)
      })
      const { records, ownTiles } = departureRecords(
        read.slice(tiles.length).flat(),
        drawn,
        keep,
      )
      await finish()
      const tooManyRead = stepBudgetError(
        rows,
        stepBudget,
        READ_STEPS_PER_KEPT,
        region,
        tile,
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
          drawn: bounds,
          attach: joinPieces(attachRows),
        })
      cut = cutRows()
      // the tiles holding the own steps of records the read lacks, each
      // fetched once for the records wanting it
      const fetches = new Map<WalkUnit, Set<string>>()
      for (const { departure } of cut.unattached) {
        const { name } = departure
        for (const unit of attachTiles(units, ownTiles.get(departure) ?? [])) {
          if (!decodedFor(unit, name)) {
            fetches.set(unit, (fetches.get(unit) ?? new Set()).add(name))
          }
        }
      }
      const wanted = [...fetches].slice(0, attachFetches)
      if (wanted.length > 0) {
        const fetched = await this.file.read(
          wanted.map(([unit]) => unit),
          reads.signal,
        )
        wanted.forEach(([unit, names], k) => {
          decode(unit, fetched[k]!, attachRows, name => names.has(name))
        })
        await finish()
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
