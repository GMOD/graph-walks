# @jbrowse/graph-walks

Cuts a window out of a pangenome graph stored as
[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) files, by range requests.

- **In:** a `gfa-to-tabix --walks` file, `<prefix>.<sample>.walks.bed.gz` with
  its `.tbi`, or a segment index, `<prefix>.segs.bed.gz` and
  `<prefix>.links.bed.gz`, and a region on the reference
- **Out:** from a walk file, the reference over the region, every haplotype's
  walk through it, and the nodes and links they visit, as typed arrays; from a
  segment index, the segments and links under the region as GFA text
- Reads walk format 4, which gfa-to-tabix writes from 1.2. The reader refuses
  any other format with a `WalkFormatError` asking for a rebuild
- Has no runtime dependencies: it reads through any object with
  [@gmod/tabix](https://github.com/GMOD/tabix-js)'s `TabixIndexedFile` methods

The JBrowse graph plugin, BandageJS, ggbandage and the sequenceTubeMap fork read
through it.

## Install

```bash
npm install @jbrowse/graph-walks @gmod/tabix
```

## Usage

```ts
import { TabixIndexedFile } from '@gmod/tabix'
import { WalkReader, walksKept } from '@jbrowse/graph-walks'
import { RemoteFile } from 'generic-filehandle2'

const url = 'https://example.org/hprc-v2.1-mc-grch38.GRCh38.walks.bed.gz'
const reader = new WalkReader(
  new TabixIndexedFile({
    filehandle: new RemoteFile(url),
    tbiFilehandle: new RemoteFile(`${url}.tbi`),
  }),
)

const refName = 'GRCh38#0#chr6'
const { graph, kept, fragments } = await reader.cut(
  refName,
  { start: 160_614_798, end: 160_647_758 },
  { keep: walksKept(['HG00097', 'HG00133'], refName) },
)
const tables = graph.tables(kept, fragments)
```

- `tables` lists the nodes with their lengths and reference positions, the links
  between them, and each walk's steps. `graphTablesGFA` in
  [@jbrowse/bandage-core](https://github.com/GMOD/bandage-core) writes it as GFA
- `reader.header()` lists the reference sample and every haplotype with rows,
  read from the walk file's first block
- `reader.referenceSequences()` lists the sequences a region can be on, in full
  PanSN: `GRCh38#0#chr6`

## A region by contig

```ts
const { tables, refName, sample } = await reader.cutRegion(
  { refName: 'chr6', start: 160_614_798, end: 160_647_758 },
  { haplotypes: ['HG00097', 'HG00133'], nodeLimit: 100_000 },
)
```

- `cutRegion` resolves `chr6` to the reference sample's sequence,
  `GRCh38#0#chr6`, and keeps the haplotypes named, every one when none are
- `referenceSample` names another sample than the header's first reference. A
  sequence the file lacks throws, naming the ones it holds
- The cut options below apply too, with `haplotypes` in place of `keep`

## Cut options

| Option       | Effect                                                                    |
| ------------ | ------------------------------------------------------------------------- |
| `keep`       | Which walks to decode, by `sample#haplotype#contig`; every one by default |
| `context`    | bp of reference kept past each end of the region; default 1000            |
| `stepBudget` | The most steps a cut keeps; no limit by default                           |
| `byteBudget` | The most compressed bytes a cut fetches; no limit by default              |
| `nodeLimit`  | The most nodes a cut keeps; no limit by default                           |
| `contained`  | `true` keeps each walk within the window; see below                       |
| `links`      | `false` lists no links, for a caller that draws no edges                  |
| `sequences`  | Keep each node's bases in `graph.sequences`, from a `--sequences` file    |
| `signal`     | An `AbortSignal` that stops the reads                                     |

- A cut over a budget or its node limit throws a `RegionTooLargeError`. The
  error's `fitsBp` holds about how wide a window would fit, and JBrowse shows it
  as a zoom-in notice
- The cut checks the node limit as it decodes each chunk, against the nodes it
  keeps whatever else it reads, and stops reading at the first chunk over it
- `WALK_BUDGETS` holds the step and byte budgets the JBrowse plugin and
  BandageJS cut by

## Fragments

Each of `fragments` is one walk's stretch through the cut: `name`, `fragStart`
(the W line's start), `hapStart`/`hapEnd`, and its steps. `before` and `after`
say what lies past each end. First come `trimmed` bp, which the cut read but did
not keep, then:

| `kind`      | The walk                                                                                                   |
| ----------- | ---------------------------------------------------------------------------------------------------------- |
| `elsewhere` | goes on in `chunk`, which the cut did not read; `bp` on, it comes back, when the read holds its next piece |
| `break`     | reaches its record's end, and the contig's next record starts `bp` on, so the graph lacks those bases      |
| `end`       | reaches the end of the contig's records in the graph                                                       |

- `graph.tables(kept, fragments)` keeps the stops as `before`/`after` columns:
  `kinds` indexes `kindNames`, `bp` is NaN where a stop states none, `trimmed`
  is 0 where it states none, and `chunkRefs` indexes `chunkRefNames`, -1 for
  none, beside `chunkStarts`
- A cut reads on to at most four chunks past what it read, so a walk that goes
  farther away and comes back stays two fragments, each `elsewhere`

## Contained cuts

- By default a cut follows each haplotype that leaves the window until it
  rejoins the reference, and keeps the reference out to there
- `contained: true` keeps the reference over the window and its context alone,
  and splits a walk wherever it steps onto the reference off that stretch. A
  detour that comes back to the window stays whole, and the cut keeps `context`
  bp of any other
- Each split leaves two fragments, and each counts the bp of the walk the cut
  read past it as `trimmed`
- A contained cut reads on to the chunk a cut-off walk goes on in only when that
  chunk is under the window and its context
- At AMY1 haplotypes go from a 4 kb window to the other copies of the gene and
  back. The default cut keeps 9,928 nodes there, and the contained cut 182

## Naming the window

- `reader.sequenceFor(sample, refName)` gives the walk file's sequence for a
  contig on a sample: `chr6` on `GRCh38` is `GRCh38#0#chr6`, and a pggb file's
  `chm13#chr1` is `chr1` on `chm13`
- `walksKept(haplotypes, refName)` is the `keep` filter for a cut on `refName`:
  its own walk and those of the haplotypes named
- `refNameLookup` and `resolveRefName` resolve a contig the same way against any
  tabix file's sequence names

## Segment indexes

```ts
import { rgfaCut } from '@jbrowse/graph-walks'

const gfa = await rgfaCut(
  { segments: segsFile, links: linksFile },
  { refName: 'GRCh38#0#chr6', start: 160_614_798, end: 160_647_758 },
  { hops: 1 },
)
```

- `segsFile` and `linksFile` are `TabixIndexedFile`s, or any object with its
  `getLines`
- gfa-to-tabix's default, anchored layout files every allele under the reference
  interval its bubble hangs from, so one read of each file holds the whole cut
- Its contig layout files each segment under the segment's own sequence. Each of
  `hops` rounds follows links out from the window, and a last read closes every
  allele reached
- `rgfaCut` writes each segment as `*` with LN/SN/SO/SR tags and the tag
  column's tags, sorted, so one region always gives the same text
- `parseSegmentRow`, `parseLinkLine` and `segmentSamples` read single rows

## How a cut works

- A chunk holds the reference's path through it once, each way a haplotype
  leaves that path once, and every other path's piece as the list of ways it
  takes; gfa-to-tabix's
  [walks-format.md](https://github.com/GMOD/gfa-to-tabix/blob/main/docs/walks-format.md)
  describes the rows
- A cut reads every chunk under the region, and the chunks before it as far as
  the header's `maxnode` reaches, since a node is filed under the chunk holding
  its start, and decodes each chunk once
- A haplotype that leaves the reference near the region is followed to where it
  rejoins, and the reference is cut out as far as those rejoins reach
- A walk the read cuts off goes on in the chunk its row names (`pv:`/`nx:`
  tags); the cut reads on to that chunk, in up to three rounds

## License

MIT
