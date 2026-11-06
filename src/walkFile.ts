import { unzip } from '@gmod/bgzf-filehandle'

import type { GenericFilehandle } from 'generic-filehandle2'

/**
 * Walk format 7 (gfa-to-tabix's ADR 0010): `<prefix>.<sample>.walks.gz`, a
 * series of gzip units, and its directory `<prefix>.<sample>.walks.gz.idx`.
 * The directory's first line is the data's `#walks` line plus
 * `bytes:i:<data length>`, then `#header 0 <length>` and one
 * `seq <offset> <length>` line per reference sequence, giving the unit that
 * lists the sequence's units: `#entry seq start end kind name offset length`.
 */
export const WALK_FORMAT = 7

// What the reader needs of a file: generic-filehandle2's handles are ones,
// and `stat` gives the data's length, from Content-Range on a RemoteFile's
// first read. Where it gives 0, as a RemoteFile does when CORS hides
// Content-Range, the reader skips checking the directory's `bytes`
export type WalkFilehandle = Pick<
  GenericFilehandle,
  'read' | 'readFile' | 'stat' | 'source'
>

// A walk file in a format other than 7, which the reader refuses
export class WalkFormatError extends Error {
  override name = 'WalkFormatError'
}

/**
 * One unit of a sequence's table: `start end` is a tile's true extent, a
 * departure block's span, or a detail record's pieces' span on the
 * reference; `name` a tile's nominal start, a detail's `tile:detail`, or `.`
 * for a block; `offset length` its compressed bytes in the data file
 */
export interface WalkUnit {
  refName: string
  start: number
  end: number
  kind: string
  name: string
  offset: number
  length: number
}

interface Directory {
  walks: string[]
  bytes: number | undefined
  header: { offset: number; length: number }
  tables: Map<string, { offset: number; length: number }>
}

const decoder = new TextDecoder()

function lines(bytes: Uint8Array) {
  const text = decoder.decode(bytes)
  return text.split('\n').filter(line => line !== '')
}

function integerTag(fields: string[], tag: string) {
  const field = fields.find(f => f.startsWith(`${tag}:i:`))
  return field === undefined ? undefined : +field.slice(tag.length + 3)
}

function formatError(format: number) {
  return new WalkFormatError(
    format > WALK_FORMAT
      ? `This walk file is format ${format}; this version of @jbrowse/graph-walks reads format ${WALK_FORMAT}, so update the viewer, or @jbrowse/graph-walks, to one that reads format ${format}`
      : `This walk file is format ${format}; @jbrowse/graph-walks reads format ${WALK_FORMAT}, so rebuild it with gfa-to-tabix 4.0 or later`,
  )
}

function parseDirectory(text: string): Directory {
  const [first = '', ...rest] = text.split('\n').filter(line => line !== '')
  const fields = first.split('\t')
  if (fields[0] !== '#walks') {
    throw new WalkFormatError(
      `This is not a walk file's directory: its first line starts ${JSON.stringify(first.slice(0, 40))}, not #walks`,
    )
  }
  const format = integerTag(fields, 'format') ?? 1
  if (format !== WALK_FORMAT) {
    throw formatError(format)
  }
  let header: Directory['header'] | undefined
  const tables = new Map<string, { offset: number; length: number }>()
  for (const line of rest) {
    const [name = '', offset, length] = line.split('\t')
    const range = { offset: Number(offset), length: Number(length) }
    if (
      !Number.isSafeInteger(range.offset) ||
      !Number.isSafeInteger(range.length)
    ) {
      throw new WalkFormatError(
        `A walk file directory line holds no byte range: ${line}`,
      )
    }
    if (name === '#header') {
      header = range
    } else {
      tables.set(name, range)
    }
  }
  if (!header) {
    throw new WalkFormatError('The walk file directory has no #header line')
  }
  return {
    walks: fields.filter(f => !f.startsWith('bytes:i:')),
    bytes: integerTag(fields, 'bytes'),
    header,
    tables,
  }
}

// A read shared by every caller, tried again by the next caller once it fails
export function once<T>(read: () => Promise<T>) {
  let held: Promise<T> | undefined
  return () => {
    if (!held) {
      const started = read()
      started.catch(() => {
        if (held === started) {
          held = undefined
        }
      })
      held = started
    }
    return held
  }
}

// Settles like `promise` unless `signal` aborts first
export function raceAbort<T>(promise: Promise<T>, signal?: AbortSignal) {
  return signal
    ? new Promise<T>((resolve, reject) => {
        signal.throwIfAborted()
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
    : promise
}

const FORMAT_6_NAME = /\.walks\.bed\.gz(?=$|[?#])/

/**
 * A format 7 walk file read by byte ranges: the directory and the header
 * once, each sequence's table once, and units as asked, the ranges of units
 * that touch merged into one single-range read
 */
export class WalkFile {
  private readonly data: WalkFilehandle
  private readonly directory: WalkFilehandle
  private readonly tableReads = new Map<string, () => Promise<WalkUnit[]>>()

  constructor(data: WalkFilehandle, directory: WalkFilehandle) {
    this.data = data
    this.directory = directory
  }

  private readonly readDirectory = once(async () => {
    try {
      return parseDirectory(await this.directory.readFile('utf8'))
    } catch (error) {
      throw !(error instanceof WalkFormatError) &&
        FORMAT_6_NAME.test(this.data.source ?? '')
        ? new WalkFormatError(
            `${this.data.source} names a walk file of format 6 or earlier, which has no .idx directory; @jbrowse/graph-walks reads format ${WALK_FORMAT}, so rebuild it with gfa-to-tabix 4.0 or later`,
            { cause: error },
          )
        : error
    }
  })

  // The header unit's lines, once the directory and the data agree
  readonly header = once(async () => {
    const directory = await this.readDirectory()
    const { offset, length } = directory.header
    const header = lines(await unzip(await this.data.read(length, offset)))
    const { size } = await this.data.stat()
    if (header[0] !== directory.walks.join('\t')) {
      throw new WalkFormatError(
        `The walk file's directory does not match its data: the directory states ${directory.walks.join(' ')}, the data ${header[0] ?? 'no #walks line'}`,
      )
    } else if (size > 0 && directory.bytes !== size) {
      throw new WalkFormatError(
        `The walk file's directory does not match its data: the directory states ${directory.bytes ?? 'no'} bytes, the data holds ${size}`,
      )
    }
    return header
  })

  // The reference sequences, in the directory's order
  async sequences() {
    return [...(await this.readDirectory()).tables.keys()]
  }

  // A sequence's units in file order; none for a sequence the file lacks
  async units(refName: string) {
    let read = this.tableReads.get(refName)
    if (!read) {
      read = once(async () => {
        const range = (await this.readDirectory()).tables.get(refName)
        return range ? this.table(refName, range) : []
      })
      this.tableReads.set(refName, read)
    }
    return read()
  }

  private async table(
    refName: string,
    { offset, length }: { offset: number; length: number },
  ) {
    const [, ...entries] = lines(
      await unzip(await this.data.read(length, offset)),
    )
    return entries.map(line => {
      const [, seq = '', start, end, kind = '', name = '', at, size] =
        line.split('\t')
      const unit = {
        refName: seq,
        start: Number(start),
        end: Number(end),
        kind,
        name,
        offset: Number(at),
        length: Number(size),
      }
      if (
        seq !== refName ||
        ![unit.start, unit.end, unit.offset, unit.length].every(
          Number.isSafeInteger,
        )
      ) {
        throw new WalkFormatError(`${refName}'s unit table holds ${line}`)
      }
      return unit
    })
  }

  /**
   * The lines of each unit, without its `#unit` line, in the order given.
   * Units whose byte ranges touch are fetched by one read, and each is
   * gunzipped on its own.
   */
  async read(units: WalkUnit[], signal?: AbortSignal) {
    const byOffset = [...new Map(units.map(u => [u.offset, u])).values()]
      .filter(u => u.length > 0)
      .sort((a, b) => a.offset - b.offset)
    const runs: WalkUnit[][] = []
    for (const unit of byOffset) {
      const run = runs.at(-1)
      const last = run?.at(-1)
      if (run && last && last.offset + last.length === unit.offset) {
        run.push(unit)
      } else {
        runs.push([unit])
      }
    }
    const held = new Map<number, string[]>()
    await Promise.all(
      runs.map(async run => {
        const first = run[0]!
        const last = run.at(-1)!
        const length = last.offset + last.length - first.offset
        const bytes = await this.data.read(length, first.offset, { signal })
        for (const unit of run) {
          const at = unit.offset - first.offset
          const [head, ...rest] = lines(
            await unzip(bytes.subarray(at, at + unit.length)),
          )
          if (!head?.startsWith('#unit\t')) {
            throw new WalkFormatError(
              `No #unit line at byte ${unit.offset} of the walk file, where ${unit.refName}'s table puts a ${unit.kind} unit`,
            )
          }
          held.set(unit.offset, rest)
        }
      }),
    )
    return units.map(unit => held.get(unit.offset) ?? [])
  }
}
