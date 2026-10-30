# @jbrowse/graph-walks

@jbrowse/graph-walks reads a window of a pangenome graph from
[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) files by range requests.

- **In:** a `gfa-to-tabix --walks` file, `<prefix>.<sample>.walks.bed.gz` with
  its `.tbi`, or a segment index, `<prefix>.segs.bed.gz` and
  `<prefix>.links.bed.gz` with theirs or format 2's `<prefix>.graph.bed.gz`,
  plus a reference region such as `chr6:160,614,798-160,647,758`
- **Out:** from a walk file, the reference over the region, every haplotype's
  walk through it, and the nodes and links they visit, as typed arrays; from a
  segment index, the segments and links under the region as GFA text
- The reader reads walk format 5, which gfa-to-tabix 1.4 and later write, and
  refuses any other format with a `WalkFormatError` asking for a rebuild
- The package has no runtime dependencies, and reads through any object with
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

const url =
  'https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38-f5s.GRCh38.walks.bed.gz'
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
- A node the reference passes more than once, as in a collapsed repeat, takes
  the position of its first pass over the window
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

| Option        | Effect                                                                              |
| ------------- | ----------------------------------------------------------------------------------- |
| `keep`        | Which walks to decode, by `sample#haplotype#contig`; every one by default           |
| `context`     | bp of reference kept past each end of the region; default 1000                      |
| `stepBudget`  | The most steps a cut keeps; no limit by default                                     |
| `byteBudget`  | The most compressed bytes a cut fetches; no limit by default                        |
| `nodeLimit`   | The most nodes a cut keeps; no limit by default                                     |
| `contained`   | `true` keeps each walk within the window; see below                                 |
| `links`       | `false` lists no links, for a caller that draws no edges                            |
| `sequences`   | Keep each node's bases in `graph.sequences`, from a `--sequences` file              |
| `signal`      | An `AbortSignal` that stops the reads                                               |
| `expand`      | Which abridged walks to read in full, by `sample#haplotype#contig`; none by default |
| `budgetNames` | The words a `RegionTooLargeError` message uses for the step and byte budgets        |

- A cut over a budget or its node limit throws a `RegionTooLargeError`. The
  error's `fitsBp` holds about how wide a window would fit, and JBrowse shows it
  as a zoom-in notice
- Over the byte budget, the error's `bytes` is what the window would fetch, so a
  viewer with a fixed window can offer to load it anyway with `byteBudget`
  raised to `bytes`. `zoomHelps` is false where no narrower window fits
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

## Abridged walks

- In a chunk over its budget, gfa-to-tabix moves the costliest paths' rows to
  detail records and leaves one `abridged` row per path. A satellite array in a
  pggb graph is the usual case
- The cut makes a separate fragment for each such row, with an `abridged` field:
  its `steps` and `bp`, its `pieces`, and the `detail` record that holds them.
  The fragment's two steps are the handles on the reference where the walk
  leaves and rejoins, and `hapStart`/`hapEnd` span its `bp`
- `graph.tables(kept, fragments)` keeps them as `abridgedSteps` and `abridgedBp`
  columns, NaN for every other walk
- With `expand` naming a walk, the cut reads its detail records and cuts it in
  full, as a file without abridging would
- `referenceSequences()` leaves out the detail records' sequences, `seq~cs`

## Haplotypes that jump over the window

- A haplotype whose deletion spans the window has no walk row in its chunks.
  gfa-to-tabix 1.5 writes a `skip` row there instead, naming the path and the
  stretch of the reference it leaves out
- `cut` returns them as `skips`: one per jump, `name`, `fragStart`, `piece`,
  `refStart` and `refEnd`, for each path that has no fragment in the cut and
  that `keep` keeps
- The writer leaves out a skip row in a chunk where the path has a piece, so the
  cut also names a jump between two reference steps it read on either side of
  the window, reading on to the far side's chunk when only the near side's is
  read, and names it even when the path has a fragment in the cut
- `graph.tables(kept, fragments, skips)` keeps them as a `skips` table of
  `names`, `starts` and `ends`, which `cutRegion` fills in
- Skip rows hold no steps or nodes, so only the byte budget counts them, with
  the rest of their chunk
- Where a haplotype leaves the reference before the window and rejoins after it,
  the cut keeps its steps off the reference between as a fragment, which counts
  the walk read past either end as `trimmed`

## Moved components

- gfa-to-tabix 1.5's `--max-component` leaves a component over the limit on its
  own coordinates, so an anchored cut reaches it only through the links into it
- `rgfaCut` tags each segment of such a component that it holds `mv:i:1`
  (`MOVED_TAG`), for a viewer to draw as a stub

## Structural forms

```ts
import {
  describeForm,
  formTag,
  formsCut,
  structuralForms,
} from '@jbrowse/graph-walks'

const region = { refName: 'chr1', start: 109_680_000, end: 109_715_000 }
const haplotypes = await formsCut(reader, region)
const result = structuralForms(haplotypes, { window: region })
result.forms.map(f => `${f.members.length} ${describeForm(f)}`)
// ['306 del 18.4kb', '156 ref']
result.forms.map(f => formTag(f, result))
// ['del 18.4kb', 'ref']
```

- `formsCut` cuts the region uncontained and without links, reads abridged walks
  in full while the cut stays under `FORMS_STEP_BUDGET` steps, and returns each
  haplotype's departures from the reference: exact REF span and ALT length,
  inverted bp, and bp walked again by a step back along the reference
- `formsCut` with `partial: false` reads every abridged walk or none: where
  reading them all would take the cut past `stepBudget`, it throws
  `ExpansionTooLargeError`, a `RegionTooLargeError` with the step count, before
  the second cut
- `formsCut` takes `byteBudget`, and a region the index estimates past it throws
  `RegionTooLargeError` before any read. The budget covers the region's chunks,
  not the detail records an expanded cut reads, which only the step budget
  bounds
- `walkDepartures` reads the same departures off any cut's fragments and
  `graph.nodes`, or off its tables
- `structuralForms` groups the haplotypes that make the same changes at every
  informative site. Each change is a `Mark`: an insertion, deletion,
  replacement, inversion, shift, duplication or `net`, the last for an abridged
  walk left unexpanded
- Each mark in a result names its `site` in `sites` and its `state` key there. A
  site is `edge`, `shared` (no form differs there) or informative, and
  `formStateAt(result, form, site)` gives a form's state key at any of them.
  `REFERENCE_STATE`, `NO_CALL_STATE` and `RARER_STATE` are the keys for no
  change, no call and a rare change among no-calls
- `placed` counts the haplotypes in a form, so a form's share is its members
  over `placed`
- `form.onlyShifts` is true where every change outside the shared sites is a
  shift and every site has a call: the reference's length with a stretch aligned
  elsewhere, as in AMY1's H3r
- `markTag`, `formTag` and `describeForm` label marks and forms, as
  [Labels](#labels) shows
- Every result holds plain data, so a worker can post it, and regrouping at
  another floor needs no new cut
- [docs/structural-forms.md](docs/structural-forms.md) gives the grouping rules
  and the measurements behind them

| Option          | Effect                                                                                   |
| --------------- | ---------------------------------------------------------------------------------------- |
| `window`        | The region; marks must touch it, and `auto` reads its width                              |
| `floorBp`       | The smallest change a form tells apart; `auto` gives 50 bp, 1 kb or 10 kb by width       |
| `minCarriers`   | The fewest haplotypes a state or form needs; 1% of them, at least one, by default        |
| `joinBp`        | Departures of one haplotype nearer than this on the reference join; default 50           |
| `exclude`       | Samples to leave out, such as a second reference carried as a haplotype; none by default |
| `foldRareForms` | Fold states and forms under `minCarriers` into the nearest common one; on by default     |

### Labels

- `markTag(mark)` is the per-mark label and takes a `Mark` or a `SiteMark`. C4's
  common change gives `del 6.46kb->94bp`
- `formTag` names a form in about 20 characters for a row label: its largest
  change's `markTag`, then ` +N` counting the rest. Given the result, it leaves
  out the changes every form shares
- `describeForm` lists every change's `markTag`, largest first, comma-separated
- Sizes take 3 significant figures with the unit attached (`94bp`, `6.46kb`,
  `1.2Mb`), as `sizeTag` writes them
- A replacement whose smaller side is under 5% (`LOPSIDED_RATIO`) of its larger
  side stays `kind: 'replacement'` with its exact REF span and ALT length, and
  gains `lopsided: 'deletion' | 'insertion'`, so a caller can colour it as a
  loss or gain. Its label reads as one

| Change                         | Token                                   |
| ------------------------------ | --------------------------------------- |
| deletion, insertion            | `del 18.4kb`, `ins 22kb`                |
| lopsided replacement           | `del 6.46kb->94bp`, `ins 178bp->13.9kb` |
| other replacement              | `repl 32.7kb->58.7kb`                   |
| shift, inversion, duplication  | `shift 94.2kb`, `inv 12kb`, `dup 49kb`  |
| shift that changes length      | `shift 94.2kb net +2.81kb`              |
| abridged walk left unexpanded  | `net +78kb`, `net -94kb`                |
| no change, only shared changes | `ref`, `shared`                         |
| rare change, no call           | `rarer`, `no call`                      |
| more changes, in `formTag`     | `del 32.7kb +1`                         |

## Contained cuts

- By default a cut follows each haplotype that leaves the region near it until
  it rejoins the reference, and keeps the reference out as far as those rejoins
  reach
- `contained: true` keeps the reference over the window and its context alone,
  and splits a walk wherever it steps onto the reference off that stretch. A
  detour that comes back to the window stays whole, and the cut keeps `context`
  bp of any other
- Each split leaves two fragments, and each counts the bp of the walk the cut
  read past it as `trimmed`
- A contained cut reads on to the chunk a cut-off walk goes on in only when that
  chunk is under the window and its context, or holds the far side of a walk
  that crosses the window
- At AMY1 haplotypes go from a 4 kb window to the other copies of the gene and
  back. The default cut keeps 9,928 nodes there, and the contained cut 182

## Naming the window

- `reader.sequenceFor(sample, refName)` gives the walk file's sequence for a
  contig on a sample: `chr6` on `GRCh38` is `GRCh38#0#chr6`, and a pggb file's
  `chm13#chr1` is `chr1` on `chm13`
- `walksKept(haplotypes, refName)` returns the `keep` filter for a cut on
  `refName`, which keeps the reference walk on `refName` and the walks of the
  haplotypes named
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
- For gfa-to-tabix's format 2, pass `{ graph: graphFile }`, the one
  `<prefix>.graph.bed.gz`, with `getLines` and `getHeaderLines`. `rgfaCut` reads
  its layout from the header, refuses a format other than 2, and writes each
  node's `SB:H:` bits as the `SM:Z:` names format 1 states, in `#samples` order.
  A link end its row gives no carriers takes them from its node row in the same
  query, so a cut of either format gives the same segments and links
- gfa-to-tabix's default, anchored layout files every allele under the reference
  interval its bubble hangs from, so one read of each file holds the whole cut
- Its contig layout files each segment under the segment's own sequence. Each of
  `hops` rounds follows links out from the window, and a last read closes every
  allele reached
- `rgfaCut` writes each segment as `*` with LN/SN/SO/SR tags and the tag
  column's tags, sorted, so one region always gives the same text
- `parseSegmentRow`, `parseLinkLine` and `segmentSamples` read single rows;
  `parseGraphHeader`, `parseGraphRow` and `carrierNames` read format 2's

## How a cut works

- A chunk holds the reference's path through it once, each way a haplotype
  leaves that path once, and every other path's piece as the list of ways it
  takes; gfa-to-tabix's
  [walks-format.md](https://github.com/GMOD/gfa-to-tabix/blob/main/docs/walks-format.md)
  describes the rows
- A cut reads every chunk under the region, from the chunk before the region's
  first, or the chunk before a `#long` pass that overlaps the region, whichever
  is earlier. It decodes each chunk once
- A file with no `#long` lines, from gfa-to-tabix 1.2, makes the cut read
  further back when the header's `maxnode` exceeds a chunk
- A walk the read cuts off goes on in the chunk its row names (`pv:`/`nx:`
  tags). The cut reads on to that chunk, in up to three rounds and at most four
  chunks past what it read

## License

MIT
