import { panSNHaplotype, panSNSample } from './util.ts'

import type { GraphTables } from './util.ts'
import type { WalkReader, WalkReaderHeader } from './walkReader.ts'
import type { WalkFragment, WalkNode } from './walkRows.ts'

/**
 * Where a haplotype leaves the reference path R and rejoins it: R's bases
 * `refStart`–`refEnd` against `altLength` of its own. `invertedBp` counts the
 * steps along R against its orientation. A step back along R walks a tandem
 * copy again: the departure is then a point at the junction, `repeatedBp`
 * counts the bp walked again, and `copied` spans them on R. An `abridged`
 * departure stands for the pieces a writer abridged, so only its net change
 * is exact.
 */
export interface Departure {
  refStart: number
  refEnd: number
  altLength: number
  invertedBp: number
  repeatedBp: number
  copied?: [number, number]
  abridged?: true
}

/**
 * A `sample#haplotype`'s departures, by refStart. `covered` spans the R its
 * walks reach, and `open` says one of them leaves R and never rejoins in the
 * cut, so it has no call past that end.
 */
export interface HaplotypeDepartures {
  haplotype: string
  covered: [number, number][]
  departures: Departure[]
  open: boolean
}

export type DepartureSource =
  { fragments: WalkFragment[]; nodes: Map<number, WalkNode> } | GraphTables

interface Walk {
  name: string
  hapStart: number
  hapEnd: number
  steps: number
  id: (k: number) => number
  rev: (k: number) => boolean
  abridgedBp?: number
  // the walk's bp a cut read but kept out of this one, before and after it
  trimmed: [number, number]
}

interface Lookup {
  walks: Walk[]
  // R position of a node on R, NaN off it
  refStart: (id: number) => number
  length: (id: number) => number
}

function fromFragments(
  fragments: WalkFragment[],
  nodes: Map<number, WalkNode>,
  refName: string,
): Lookup {
  return {
    walks: fragments.map(f => ({
      name: f.name,
      hapStart: f.hapStart,
      hapEnd: f.hapEnd,
      steps: f.ids.length,
      id: k => f.ids[k]!,
      rev: k => f.rev[k] === 1,
      abridgedBp: f.abridged?.bp,
      trimmed: [f.before.trimmed ?? 0, f.after.trimmed ?? 0],
    })),
    refStart: id => {
      const node = nodes.get(id)
      return node?.rank === 0 && node.refName === refName ? node.start : NaN
    },
    length: id => {
      const node = nodes.get(id)
      return node ? node.end - node.start : 0
    },
  }
}

function fromTables({ nodes, walks }: GraphTables, refName: string): Lookup {
  const ref = nodes.refNames.indexOf(refName)
  return {
    walks: walks.names.map((name, p) => {
      const from = walks.offsets[p]!
      const bp = walks.abridgedBp?.[p]
      return {
        name,
        hapStart: walks.starts[p]!,
        hapEnd: walks.ends[p]!,
        steps: walks.offsets[p + 1]! - from,
        id: k => walks.steps[from + k]!,
        rev: k => walks.reversed[from + k] === 1,
        abridgedBp: bp === undefined || Number.isNaN(bp) ? undefined : bp,
        trimmed: [0, 0],
      }
    }),
    refStart: id =>
      nodes.ranks[id] === 0 && nodes.refs[id] === ref ? nodes.starts[id]! : NaN,
    length: id => nodes.lengths[id] ?? 0,
  }
}

const onR = (start: number) => !Number.isNaN(start)

/**
 * A path's walks in order, those that abut on the haplotype joined. An
 * abridged walk's pieces need not be consecutive: where the path steps onto
 * another chunk's R and back, the pieces between lie in that chunk. A walk
 * starting before the abridged pieces' bp run out is one of those, so its
 * whole run, trimmed bp included, joins the abridged walk's bp.
 */
function chains(walks: Walk[]) {
  const sorted = walks.toSorted((a, b) =>
    a.name < b.name ? -1 : a.name > b.name ? 1 : a.hapStart - b.hapStart,
  )
  const out: Walk[][] = []
  let end = NaN
  let runEnd = -Infinity
  for (const walk of sorted) {
    const chain = out.at(-1)
    const last = chain?.at(-1)
    if (chain && last?.name === walk.name) {
      if (last.abridgedBp !== undefined && walk.hapStart < end) {
        const start = walk.hapStart - walk.trimmed[0]
        const stop = walk.hapEnd + walk.trimmed[1]
        const added = Math.max(0, stop - Math.max(start, runEnd))
        runEnd = Math.max(runEnd, stop)
        end += added
        chain[chain.length - 1] = {
          ...last,
          abridgedBp: last.abridgedBp + added,
        }
        continue
      }
      if (walk.hapStart === end) {
        chain.push(walk)
        end = walk.hapEnd
        runEnd = -Infinity
        continue
      }
    }
    out.push([walk])
    end = walk.hapEnd
    runEnd = -Infinity
  }
  return out
}

/**
 * Each haplotype's departures from R, `refName`, read off a cut's fragments
 * and `graph.nodes`, or off its tables. Fragments of a path that abut on the
 * haplotype read as one walk, and one that steps along R mostly reversed
 * reads reversed. Steps off R before a walk's first forward step on R or
 * after its last leave it open there. An abridged fragment's bp, with those
 * of any pieces between its own that lie in other chunks, lie between the R
 * steps either side of it, so the departure holding it is exact in net. At
 * either end of a walk it has no R step beyond it and stops at its anchors,
 * which is exact only where the walk enters or leaves R there. The walks on
 * R's own haplotype are left out.
 */
export function walkDepartures(source: DepartureSource, refName: string) {
  const { walks, refStart, length } =
    'fragments' in source
      ? fromFragments(source.fragments, source.nodes, refName)
      : fromTables(source, refName)
  const reference = panSNHaplotype(refName) ?? refName
  const byHaplotype = new Map<string, HaplotypeDepartures>()
  for (const chain of chains(walks)) {
    const { name } = chain[0]!
    const haplotype = panSNHaplotype(name) ?? panSNSample(name)
    if (haplotype === reference || name === refName) {
      continue
    }
    let h = byHaplotype.get(haplotype)
    if (!h) {
      h = { haplotype, covered: [], departures: [], open: false }
      byHaplotype.set(haplotype, h)
    }
    readChain(chain, refStart, length, h)
  }
  const out = [...byHaplotype.values()]
  for (const h of out) {
    h.departures.sort((a, b) => a.refStart - b.refStart || a.refEnd - b.refEnd)
    h.covered = mergeSpans(h.covered)
  }
  return out
}

function readChain(
  chain: Walk[],
  refStart: (id: number) => number,
  length: (id: number) => number,
  h: HaplotypeDepartures,
) {
  let forwardBp = 0
  let reversedBp = 0
  for (const walk of chain) {
    for (let k = 0; walk.abridgedBp === undefined && k < walk.steps; k++) {
      const id = walk.id(k)
      if (onR(refStart(id))) {
        if (walk.rev(k)) {
          reversedBp += length(id)
        } else {
          forwardBp += length(id)
        }
      }
    }
  }
  const flip = reversedBp > forwardBp
  let first = NaN
  let prevEnd = NaN
  let reach = -Infinity
  let alt = 0
  let inverted = 0
  // set while an abridged fragment's bp are in `alt`: its anchors' far end
  // on R, and `alt` there
  let abridgedEnd = NaN
  let altAtAbridgedEnd = 0
  for (const walk of flip ? chain.toReversed() : chain) {
    const n = walk.steps
    if (walk.abridgedBp !== undefined) {
      const a = walk.id(0)
      const b = walk.id(n - 1)
      const sa = refStart(a)
      const sb = refStart(b)
      if (!onR(sa) || !onR(sb)) {
        h.open = true
        continue
      }
      if (Number.isNaN(prevEnd)) {
        first = Math.min(sa, sb)
        prevEnd = first
        alt = 0
      }
      alt += walk.abridgedBp
      abridgedEnd = Math.max(sa + length(a), sb + length(b))
      altAtAbridgedEnd = alt
      continue
    }
    for (let k = 0; k < n; k++) {
      const step = flip ? n - 1 - k : k
      const id = walk.id(step)
      const reversed = walk.rev(step) !== flip
      const start = refStart(id)
      const len = length(id)
      if (onR(start) && !reversed) {
        if (Number.isNaN(prevEnd)) {
          first = start
        } else if (!Number.isNaN(abridgedEnd)) {
          h.departures.push({
            refStart: prevEnd,
            refEnd: Math.max(start, prevEnd),
            altLength: alt + Math.max(0, prevEnd - start),
            invertedBp: 0,
            repeatedBp: 0,
            abridged: true,
          })
          abridgedEnd = NaN
        } else if (start !== prevEnd || alt > 0) {
          const repeated = start < prevEnd ? prevEnd - start : 0
          h.departures.push({
            refStart: prevEnd,
            refEnd: repeated > 0 ? prevEnd : start,
            altLength: alt + repeated,
            invertedBp: inverted,
            repeatedBp: repeated,
            ...(repeated > 0 ? { copied: [start, prevEnd] } : {}),
          })
        }
        prevEnd = start + len
        reach = Math.max(reach, prevEnd)
        alt = 0
        inverted = 0
      } else if (Number.isNaN(prevEnd)) {
        h.open = true
      } else {
        alt += len
        if (onR(start)) {
          inverted += len
        }
      }
    }
  }
  if (!Number.isNaN(abridgedEnd)) {
    h.departures.push({
      refStart: prevEnd,
      refEnd: Math.max(abridgedEnd, prevEnd),
      altLength: altAtAbridgedEnd,
      invertedBp: 0,
      repeatedBp: 0,
      abridged: true,
    })
    reach = Math.max(reach, abridgedEnd)
    alt -= altAtAbridgedEnd
  }
  if (Number.isNaN(first)) {
    h.open = true
  } else {
    h.covered.push([first, reach])
  }
  if (alt > 0) {
    h.open = true
  }
}

function mergeSpans(spans: [number, number][]) {
  const sorted = spans.toSorted((a, b) => a[0] - b[0])
  const merged: [number, number][] = []
  for (const [start, end] of sorted) {
    const last = merged.at(-1)
    if (last && start <= last[1]) {
      last[1] = Math.max(last[1], end)
    } else {
      merged.push([start, end])
    }
  }
  return merged
}

// The most steps a forms cut keeps, expanded abridged walks included; an
// expanded cut at the MHC keeps 9.4 M
export const FORMS_STEP_BUDGET = 16_000_000

export interface FormsCutOptions {
  // read abridged walks in full, as far as the step budget allows; true by
  // default
  expand?: boolean
  stepBudget?: number
  // the most compressed bytes either cut fetches from the region's chunks,
  // past which the cut throws RegionTooLargeError; no limit by default. Only
  // stepBudget bounds the detail records expansion reads
  byteBudget?: number
  signal?: AbortSignal
}

/**
 * Every haplotype's departures over a region, for structuralForms. The cut is
 * uncontained, so a departure crossing an edge is followed to where it
 * rejoins, and lists no links. Abridged walks are read in full, those with
 * the fewest steps first, while the cut stays under `stepBudget`; any left
 * keep their net change only. A haplotype the header names with no walk in
 * the cut is listed with nothing covered. A region the index estimates past
 * `byteBudget` throws the reader's RegionTooLargeError before any read.
 */
export async function formsCut(
  reader: Pick<WalkReader, 'header' | 'sequenceFor' | 'cut'>,
  region: { refName: string; start: number; end: number },
  opts: FormsCutOptions = {},
) {
  const {
    expand = true,
    stepBudget = FORMS_STEP_BUDGET,
    byteBudget,
    signal,
  } = opts
  const header: WalkReaderHeader = await reader.header({ signal })
  const sample = header.references[0] ?? ''
  const refName =
    (await reader.sequenceFor(sample, region.refName, { signal })) ??
    region.refName
  const cutOpts = { links: false, byteBudget, signal }
  let cut = await reader.cut(refName, region, cutOpts)
  if (expand) {
    const expanded = expandedWalks(cut.fragments, stepBudget)
    if (expanded.size > 0) {
      cut = await reader.cut(refName, region, {
        ...cutOpts,
        expand: name => expanded.has(name),
      })
    }
  }
  const haplotypes = walkDepartures(
    { fragments: cut.fragments, nodes: cut.graph.nodes },
    refName,
  )
  const seen = new Set(haplotypes.map(h => h.haplotype))
  for (const haplotype of header.haplotypes) {
    if (!seen.has(haplotype)) {
      haplotypes.push({ haplotype, covered: [], departures: [], open: false })
    }
  }
  return haplotypes
}

// The abridged walks a cut can read in full within `stepBudget`, fewest
// steps first
function expandedWalks(fragments: WalkFragment[], stepBudget: number) {
  let kept = 0
  const abridged = new Map<string, number>()
  for (const f of fragments) {
    if (f.abridged) {
      abridged.set(f.name, (abridged.get(f.name) ?? 0) + f.abridged.steps)
    } else {
      kept += f.ids.length
    }
  }
  const expanded = new Set<string>()
  for (const [name, steps] of [...abridged].sort((a, b) => a[1] - b[1])) {
    if (kept + steps > stepBudget) {
      break
    }
    kept += steps
    expanded.add(name)
  }
  return expanded
}
