import {
  addLink,
  anchoredCut,
  closingLinks,
  closingSpans,
  formatSubgraph,
  parseGraphHeader,
  parseGraphRow,
  parseLinkLine,
  parseSegmentRow,
} from './rgfaBed.ts'
import { checkRange } from './walkReader.ts'

import type {
  GraphHeader,
  RgfaLink,
  RgfaSegment,
  StableSpan,
} from './rgfaBed.ts'
import type { TabixSource } from './walkReader.ts'

// What a hop follows: alleles, never the backbone. A rank-0 segment reached
// through a link is flanking backbone, usually outside the window, so querying
// its links walks the reference away from the region, brings back more backbone,
// and closes nothing. Every segment that closes a bubble is an allele interior,
// and an allele interior is rank > 0 in both index formats (rGFA's SR tag, and
// pggb's rank 0-or-1 against the anchor path). Reached backbone still becomes a
// node and its edge still draws; it is only not walked through.
//
// This is gfabase's rule, arrived at from the other end. `gfabase sub` expands
// from a seed either to the whole connected component (`--connected`, which its
// own README calls overkill) or with `--cutpoints 1`, which takes the biconnected
// component by *stopping at cutpoints*, the segments every end-to-end walk of
// the chromosome must traverse. In a reference-anchored pangenome the cutpoints
// are the backbone, and rank is what the index already states about it, so
// `rank > 0` is that stopping rule without a connectivity index to build. The
// third tool is the same shape again: `gfatools view -R` maps the region onto its
// precomputed bubble decomposition and returns whole overlapping bubbles
// (gfa_query_by_reg), which is complete for the same reason. Bandage's own
// `--distance` is a node-step radius, but Bandage holds the entire graph in
// memory, so its scope is a drawing filter rather than a fetch strategy.
//
// Measured with this rule: the E. coli paa window converges at one hop (17
// segments at 1 and at 2), so the cut closes and stays closed. HPRC's amylase
// window keeps growing (63 at 0, 78 at 1, 92 at 2) because its alleles have
// alleles, and wall-clock is flat across all three, since a hop's queries go out
// together and the remote index dominates.
//
// The exact version, when a graph ships a bubble index beside its segments:
// gfatools' bubble rows carry the member segment ids of every bubble, so the
// bubbles overlapping the region state which segments a complete cut must
// contain. That turns the hop count into a termination condition rather than a
// setting, and MinigraphBubbleAdapter already parses those rows.
function offReference(segments: RgfaSegment[]) {
  return segments.filter(segment => segment.rank > 0)
}

// A gfa-to-tabix index, read as @gmod/tabix's TabixIndexedFile reads it:
// format 1's two files (gfa-to-tabix docs/adr/0002-two-index-files.md), or
// format 2's one (docs/adr/0007-node-link-format-2.md)
export type RgfaFiles =
  | {
      segments: Pick<TabixSource, 'getLines'>
      links: Pick<TabixSource, 'getLines'>
    }
  | { graph: Pick<TabixSource, 'getLines' | 'getHeaderLines'> }

const headers = new WeakMap<object, Promise<GraphHeader>>()

function graphHeader(
  file: Pick<TabixSource, 'getHeaderLines'>,
  signal?: AbortSignal,
) {
  let header = headers.get(file)
  if (!header) {
    header = file.getHeaderLines({ signal }).then(parseGraphHeader)
    header.catch(() => headers.delete(file))
    headers.set(file, header)
  }
  return header
}

function getLines(
  file: Pick<TabixSource, 'getLines'>,
  refName: string,
  start: number,
  end: number,
  opts: { signal?: AbortSignal; lineCallback: (line: string) => void },
) {
  checkRange(refName, start, end)
  return file.getLines(refName, start, end, opts)
}

// The segments and links filed over `span`, and whether the index is
// anchored. Format 1 reads the segment file only when asked, as a hop needs
// only links. Format 2's one query returns both kinds, and a link end its row
// gives no carriers takes the tags of its node row in the same query; an end
// past the span has none.
async function readSpan(
  files: RgfaFiles,
  span: StableSpan,
  withSegments: boolean,
  signal?: AbortSignal,
) {
  const segments: RgfaSegment[] = []
  const links: RgfaLink[] = []
  const read = (
    file: Pick<TabixSource, 'getLines'>,
    lineCallback: (line: string) => void,
  ) =>
    getLines(file, span.refName, span.start, span.end, { signal, lineCallback })
  if ('graph' in files) {
    const header = await graphHeader(files.graph, signal)
    await read(files.graph, line => {
      const row = parseGraphRow(line, header)
      if (row?.kind === 'N') {
        segments.push(row.segment)
      } else if (row?.kind === 'L') {
        links.push(row.link)
      }
    })
    const tags = new Map(segments.map(s => [s.id, s.tags]))
    for (const link of links) {
      link.sourceSegment.tags ||= tags.get(link.source) ?? ''
      link.targetSegment.tags ||= tags.get(link.target) ?? ''
    }
    return {
      segments: withSegments ? segments : [],
      links,
      anchored: header.anchored,
    }
  }
  let anchored = false
  await Promise.all([
    withSegments &&
      read(files.segments, line => {
        const row = parseSegmentRow(line)
        anchored ||= row.anchored
        segments.push(row.segment)
      }),
    read(files.links, line => {
      links.push(parseLinkLine(line))
    }),
  ])
  return { segments, links, anchored }
}

/**
 * The graph around a window of one stable sequence, as GFA text. The window's
 * own segments come from the segment file; the link file supplies the edges
 * incident to them *and* the coordinates of the segments on the other end,
 * which typically sit on a different stable sequence (a rank>0 bubble) and so
 * are not reachable by any coordinate query on this window.
 *
 * `hops` is extra rounds of link-following past the window's own segments and
 * their immediate neighbours, each costing one tabix query per off-reference
 * segment newly reached. After the last round, one more read per stretch of
 * stable sequence the cut holds closes it (closingSpans), so an allele reached
 * from both ends draws whole. The view asks for one hop; a hopless cut stays
 * as cheap as it was. `signal` goes to every one of those queries, so a cut
 * the view has replaced stops.
 */
export async function rgfaCut(
  files: RgfaFiles,
  window: StableSpan,
  { hops = 0, signal }: { hops?: number; signal?: AbortSignal } = {},
) {
  const segments = new Map<string, RgfaSegment>()
  const links = new Map<string, RgfaLink>()
  // Both layouts read the links over the window, so that read goes out
  // beside the segment read rather than after it.
  const first = await readSpan(files, window, true, signal)
  for (const segment of first.segments) {
    segments.set(segment.id, segment)
  }
  const windowLinks = first.links
  // An anchored index files every allele and every link under the reference
  // interval its bubble hangs from, so one read of each file holds the whole
  // graph under the window and there is nothing to hop to or close.
  if (first.anchored) {
    const cut = anchoredCut(segments, windowLinks, window)
    return formatSubgraph(cut.segments, cut.links)
  }
  const linksOver = async (span: StableSpan) =>
    (await readSpan(files, span, false, signal)).links
  const addLinksOver = async (span: StableSpan) =>
    (await linksOver(span)).flatMap(link => addLink(segments, links, link))
  let frontier = offReference(
    windowLinks.flatMap(link => addLink(segments, links, link)),
  )
  for (let hop = 0; hop < hops; hop++) {
    // One hop's queries are independent of each other, so they go out
    // together rather than one round-trip at a time, and the output is sorted
    // at the end regardless of the order they arrive in.
    const reached = await Promise.all(frontier.map(addLinksOver))
    frontier = offReference(reached.flat())
  }
  if (hops > 0) {
    const spans = closingSpans(
      [...segments.values()],
      new Set(frontier.map(segment => segment.id)),
      window,
    )
    const found = await Promise.all(spans.map(span => linksOver(span)))
    for (const link of closingLinks(found.flat(), segments, spans)) {
      addLink(segments, links, link)
    }
  }
  return formatSubgraph(segments, links)
}
