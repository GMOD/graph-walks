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
- The reader reads walk format 6, which gfa-to-tabix 3.0 and later write, and
  refuses an older format with a `WalkFormatError` asking for a rebuild, and a
  newer one asking for a newer reader
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
const { graph, kept, fragments, departures } = await reader.cut(
  refName,
  { start: 160_614_798, end: 160_647_758 },
  { keep: walksKept(['HG00097', 'HG00133'], refName) },
)
const tables = graph.tables(kept, fragments, departures)
```

- `tables` lists the nodes with their lengths and reference positions, the links
  between them, each walk's steps, and every departure from the reference over
  the window ([Departures](#departures)). `graphTablesGFA` in
  [@jbrowse/bandage-core](https://github.com/GMOD/bandage-core) writes the
  nodes, links and walks as GFA
- A node the reference passes more than once, as in a collapsed repeat, takes
  the position of its pass over the window, and each walk's step on it the pass
  its row names
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

| Option                | Effect                                                                              |
| --------------------- | ----------------------------------------------------------------------------------- |
| `keep`                | Which walks to decode, by `sample#haplotype#contig`; every one by default           |
| `context`             | bp of reference kept past each end of the region; default 1000                      |
| `stepBudget`          | The most steps a cut keeps, departures' own steps included; no limit by default     |
| `byteBudget`          | The most compressed bytes the walk read fetches; no limit by default                |
| `departureByteBudget` | The most compressed bytes the departure read fetches; no limit by default           |
| `attachBp`            | The most bp of own sequence a departure gets its steps attached for; default 65,536 |
| `attachFetches`       | The most chunks fetched for own steps the read lacks; default 8                     |
| `nodeLimit`           | The most nodes a cut keeps; no limit by default                                     |
| `contained`           | `true` keeps each walk within the window; see below                                 |
| `links`               | `false` lists no links, for a caller that draws no edges                            |
| `sequences`           | Keep each node's bases in `graph.sequences`, from a `--sequences` file              |
| `signal`              | An `AbortSignal` that stops the reads                                               |
| `expand`              | Which abridged walks to read in full, by `sample#haplotype#contig`; none by default |
| `budgetNames`         | The words a `RegionTooLargeError` message uses for the step and byte budgets        |

- A cut over a budget or its node limit throws a `RegionTooLargeError`. The
  error's `fitsBp` holds about how wide a window would fit, and JBrowse shows it
  as a zoom-in notice
- Over the byte budget, the error's `bytes` is what the window would fetch, so a
  viewer with a fixed window can offer to load it anyway with `byteBudget`
  raised to `bytes`. `zoomHelps` is false where no narrower window fits
- The cut checks the node limit as it decodes each chunk, against the nodes it
  keeps whatever else it reads, and stops reading at the first chunk over it
- A window past `departureByteBudget` throws a `DepartureReadError`, which is
  not a `RegionTooLargeError`: a departure's record is filed over its whole
  stretch, so zooming in does not lower that read. A collapsed repeat that
  haplotypes hop between the copies of is the usual cause
- `WALK_BUDGETS` holds the budgets and `attachBp` the JBrowse plugin and
  BandageJS cut by

## Fragments

Each of `fragments` is one walk's stretch through the cut: `name`, `fragStart`
(the W line's start), `hapStart`/`hapEnd`, and its steps, from the walk's first
step on the drawn stretch of the reference, the window and its context, to its
last. `trimmedBefore` and `trimmedAfter` count the bp of the walk the cut read
but kept out past either end.

- A cut reads the window's chunks and the chunk before, and nothing else: where
  a walk goes past the drawn stretch and comes back, its fragment ends at its
  last step on the stretch, and its departure says where it rejoins
- A walk's fragment therefore does not depend on which other walks the cut keeps
- `graph.tables(kept, fragments, departures)` keeps `trimmedBefore` and
  `trimmedAfter` as columns of `walks`

## Departures

A departure is where a walk's consecutive steps on the reference are not
adjacent there: it leaves the reference at one breakpoint and rejoins it at
another, with `ownBp` of its own sequence between. `cut` returns every departure
over the drawn stretch, each once, from two reads: the walk read gives those
inside one piece of a walk, which the cut computes from the piece's steps, and a
read of the sequence's `~departure` sibling gives the rest as the file's
records, which gfa-to-tabix writes for every departure between two pieces.

```ts
import { departureKind } from '@jbrowse/graph-walks'

for (const d of departures) {
  // HG00097#1#chr6 leaves 160,620,000 and rejoins 160,674,000 (forward), 120,000 bp of its own
  console.log(
    d.name,
    d.leavePos,
    d.rejoinPos,
    departureKind(d),
    d.ownBp,
    d.steps?.ids.length,
  )
}
```

- `leavePos` and `rejoinPos` are breakpoints on the reference, where the walk's
  sequence stops following it and starts again; `leaveReversed` and
  `rejoinReversed` say the step runs against the reference's orientation at its
  pass
- `hapStart` is where the own sequence starts on the walk's contig and `ownBp`
  how far it runs, 0 for a direct jump; `missingBp` is the part of it between
  the contig's W lines, which the graph lacks
- `departureKind` derives what the departure does from its ends: `switch`, the
  strands differ; `back`, the walk jumps back over reference it walks again;
  `adjacent`, sequence of its own between adjacent steps; `forward`, a jump
  forward over reference it leaves out. The tools make no further call
- A departure of at most `attachBp` of own sequence carries its own steps in
  `steps`, from the pieces the cut read; for a record whose leave piece lies
  outside the read, the cut fetches that chunk, at most `attachFetches` chunks a
  cut. `fromSteps` is true for a departure computed from a piece's steps
- Those steps' nodes are kept, so the cut draws a short insertion or replacement
  at the window's edge whole, and counts them against `stepBudget`
- `graph.tables(kept, fragments, departures)` keeps them as the `departures`
  table: `names`, `hapStarts`, `ownBp`, `missingBp`, `leavePos`,
  `leaveReversed`, `rejoinPos`, `rejoinReversed`, `fromSteps`, and each one's
  own steps as node indices at `offsets[k]` to `offsets[k + 1]` of `steps`, with
  `reversed`

## Abridged walks

- In a chunk over its budget, gfa-to-tabix moves the costliest paths' rows to
  detail records and leaves one `abridged` row per path. A satellite array in a
  pggb graph is the usual case
- The cut makes a separate fragment for each such row, with an `abridged` field:
  its `steps` and `bp`, its `pieces`, and the `detail` record that holds them.
  The fragment's two steps are the handles on the reference where the walk
  leaves and rejoins, and `hapStart`/`hapEnd` span its `bp`
- `graph.tables(kept, fragments, departures)` keeps them as `abridgedSteps` and
  `abridgedBp` columns, NaN for every other walk
- With `expand` naming a walk, the cut reads its detail records and cuts it in
  full, as a file without abridging would
- `referenceSequences()` leaves out the detail records' sequences, `seq~cs`, and
  the departure siblings, `seq~departure`

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
- `formsCut` reads each haplotype's departures off the cut's `departures`, every
  one over the region, and its fragments' ends, so a haplotype that leaves the
  region and rejoins the reference far away, or jumps over it whole, is a
  departure with no step on it; `walkDepartures` does the same for any cut
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
