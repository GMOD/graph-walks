import { crc32, deflateRawSync } from 'node:zlib'

import { BlobFile } from 'generic-filehandle2'

import type { FilehandleOptions } from 'generic-filehandle2'

import type { GraphTables } from './util.ts'
import type { WalkUnit } from './walkFile.ts'

// The GFA a cut's tables stand for, as @jbrowse/bandage-core writes it
export function graphTablesGFA({ nodes, links, walks }: GraphTables) {
  const { names } = nodes
  const lines = ['H\tVN:Z:1.1']
  for (let i = 0; i < nodes.lengths.length; i++) {
    const ref = nodes.refs[i]!
    const stable =
      ref < 0
        ? ''
        : `\tSN:Z:${nodes.refNames[ref]}\tSO:i:${nodes.starts[i]}\tSR:i:${nodes.ranks[i]}`
    lines.push(`S\t${names[i]}\t*\tLN:i:${nodes.lengths[i]}${stable}`)
  }
  for (let i = 0; i < links.from.length; i++) {
    const strands = links.strands[i]!
    lines.push(
      `L\t${names[links.from[i]!]}\t${strands & 1 ? '-' : '+'}\t${names[links.to[i]!]}\t${strands & 2 ? '-' : '+'}\t0M`,
    )
  }
  for (let p = 0; p < walks.names.length; p++) {
    const [sample, haplotype, ...contig] = walks.names[p]!.split('#')
    let body = ''
    for (let s = walks.offsets[p]!; s < walks.offsets[p + 1]!; s++) {
      body += (walks.reversed[s] ? '<' : '>') + names[walks.steps[s]!]
    }
    lines.push(
      `W\t${sample}\t${haplotype}\t${contig.join('#')}\t${walks.starts[p]}\t${walks.ends[p]}\t${body}`,
    )
  }
  return lines.join('\n')
}

const BGZF_BLOCK = 0xff00
const hex = (text: string) =>
  Uint8Array.from(text.match(/../g) ?? [], byte => parseInt(byte, 16))
const BGZF_HEADER = hex('1f8b08040000000000ff060042430200')
const BGZF_EOF = hex('1f8b08040000000000ff0600424302001b0003000000000000000000')

// `text` as BGZF blocks, as gfa-to-tabix writes a unit
export function bgzf(text: string) {
  const data = new TextEncoder().encode(text)
  const blocks: Uint8Array<ArrayBuffer>[] = []
  for (let at = 0; at < data.length; at += BGZF_BLOCK) {
    const piece = data.subarray(at, at + BGZF_BLOCK)
    const deflated = deflateRawSync(piece)
    const block = new Uint8Array(26 + deflated.length)
    block.set(BGZF_HEADER)
    block.set(deflated, 18)
    const view = new DataView(block.buffer)
    view.setUint16(16, block.length - 1, true)
    view.setUint32(18 + deflated.length, crc32(piece), true)
    view.setUint32(22 + deflated.length, piece.length, true)
    blocks.push(block)
  }
  return blocks
}

const lengthOf = (blocks: Uint8Array[]) =>
  blocks.reduce((sum, block) => sum + block.length, 0)

// A file whose bytes past `from` lie `gap` further on than the blob holds them
class GappedFile extends BlobFile {
  private readonly from: number
  private readonly gap: number

  constructor(blob: Blob, from: number, gap: number) {
    super(blob)
    this.from = from
    this.gap = gap
  }

  override read(length: number, position = 0, opts?: FilehandleOptions) {
    return super.read(
      length,
      position < this.from ? position : position - this.gap,
      opts,
    )
  }

  override async stat() {
    const { size } = await super.stat()
    return { size: size + this.gap }
  }
}

const spanOf = (cols: string[]) =>
  (cols[9] ?? '').split(',').reduce((end, n) => end + Number(n), +cols[5]!)

/**
 * A format 7 walk file in memory from `lines`, each space-separated: the
 * header's `#` lines; tile rows as `seq tile kind ...`, a tile's rows from
 * `kind` on under its nominal start; and departure lines as the file holds
 * them, `seq start end departure ...`, each filed in every block of the
 * tile size it overlaps. A tile's extent is its ref rows' span, and each of
 * its lines states that extent; every inline tile has ref rows. `reversed`
 * writes a sequence's units last first, and `gap` puts that many bytes the
 * file never reads between the header and the first unit.
 */
export function walkFileOf(
  lines: string[],
  { reversed = false, gap = 0 } = {},
) {
  const rows = lines.map(line => line.split(' '))
  const header = rows.filter(r => r[0]!.startsWith('#')).map(r => r.join('\t'))
  const tile = Number(/tile:i:(\d+)/.exec(header[0] ?? '')?.[1] ?? 65536)
  const body = rows.filter(r => !r[0]!.startsWith('#'))
  const sequences = [...new Set(body.map(r => r[0]!))]
  const parts = bgzf(`${header.join('\n')}\n`)
  const headerLength = lengthOf(parts)
  let offset = headerLength + gap
  const directory = [`#header\t0\t${headerLength}`]
  const units: WalkUnit[] = []
  const write = (text: string) => {
    const blocks = bgzf(text)
    const length = lengthOf(blocks)
    for (const block of blocks) {
      parts.push(block)
    }
    offset += length
    return { offset: offset - length, length }
  }
  for (const seq of sequences) {
    const items: { start: number; order: number; unit: () => WalkUnit }[] = []
    const tiles = new Map<number, string[][]>()
    const blocks = new Map<number, string[]>()
    for (const r of body.filter(r => r[0] === seq)) {
      if (Number.isNaN(Number(r[2]))) {
        tiles.set(+r[1]!, [...(tiles.get(+r[1]!) ?? []), r])
      } else {
        for (
          let at = Math.floor(+r[1]! / tile) * tile;
          at < +r[2]!;
          at += tile
        ) {
          blocks.set(at, [...(blocks.get(at) ?? []), r.join('\t')])
        }
      }
    }
    for (const [start, tileRows] of tiles) {
      const refs = tileRows.filter(r => r[2] === 'ref')
      const lo = refs.reduce((at, r) => Math.min(at, +r[5]!), Infinity)
      const hi = refs.reduce((at, r) => Math.max(at, spanOf(r)), -Infinity)
      items.push({
        start,
        order: 0,
        unit: () => {
          const text = [
            ['#unit', seq, lo, hi, 'tile', start],
            ...tileRows.map(r => [seq, lo, hi, ...r.slice(2)]),
          ]
          const at = write(`${text.map(c => c.join('\t')).join('\n')}\n`)
          const unit = { refName: seq, start: lo, end: hi, kind: 'tile' }
          return { ...unit, name: String(start), ...at }
        },
      })
    }
    for (const [start, blockLines] of blocks) {
      items.push({
        start,
        order: 1,
        unit: () => {
          const head = ['#unit', seq, start, start + tile, 'departures', tile]
          const at = write(`${[head.join('\t'), ...blockLines].join('\n')}\n`)
          const unit = { refName: seq, start, end: start + tile }
          return { ...unit, kind: 'departures', name: String(tile), ...at }
        },
      })
    }
    items.sort((a, b) =>
      reversed
        ? b.start - a.start || b.order - a.order
        : a.start - b.start || a.order - b.order,
    )
    const held = items.map(item => item.unit())
    for (const unit of held) {
      units.push(unit)
    }
    const end = held.reduce((at, u) => Math.max(at, u.end), 0)
    const table = write(
      `${[
        ['#unit', seq, 0, end, 'table', '.'],
        ...held.map(u => [
          '#entry',
          seq,
          u.start,
          u.end,
          u.kind,
          u.name,
          u.offset,
          u.length,
        ]),
      ]
        .map(c => c.join('\t'))
        .join('\n')}\n`,
    )
    directory.push(`${seq}\t${table.offset}\t${table.length}`)
  }
  parts.push(BGZF_EOF)
  const index = `${header[0]}\tbytes:i:${lengthOf(parts) + gap}\n${directory.join('\n')}\n`
  return {
    data: new GappedFile(new Blob(parts), headerLength, gap),
    directory: new BlobFile(new Blob([index])),
    units,
  }
}
