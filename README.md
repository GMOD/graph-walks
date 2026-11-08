# @jbrowse/graph-walks

@jbrowse/graph-walks reads a window of a pangenome graph from
[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) files by range requests.

- **In:** a `gfa-to-tabix --walks` file, `<prefix>.<sample>.walks.gz` with its
  directory `<prefix>.<sample>.walks.gz.idx`, or a segment index,
  `<prefix>.segs.bed.gz` and `<prefix>.links.bed.gz` with theirs or format 2's
  `<prefix>.graph.bed.gz`, plus a reference region such as
  `chr6:160,614,798-160,647,758`
- **Out:** from a walk file, the reference over the region, every haplotype's
  walk through it, and the nodes and links they visit, as typed arrays; from a
  segment index, the segments and links under the region as GFA text
- The reader reads walk format 7, which gfa-to-tabix 4.0 and later write, and
  refuses an older format with a `WalkFormatError` asking for a rebuild, and a
  newer one asking for a newer reader
- It reads a walk file through
  [generic-filehandle2](https://github.com/GMOD/generic-filehandle2)'s handles
  and gunzips with
  [@gmod/bgzf-filehandle](https://github.com/GMOD/bgzf-filehandle); a segment
  index through any object with
  [@gmod/tabix](https://github.com/GMOD/tabix-js)'s `TabixIndexedFile` methods

The JBrowse graph plugin, BandageJS, ggbandage and the sequenceTubeMap fork read
through it.

## Install

```bash
npm install @jbrowse/graph-walks
```

A segment index also needs `@gmod/tabix`.

## Usage

```ts
import { WalkReader, walksKept } from '@jbrowse/graph-walks'

const reader = WalkReader.open(
  'https://example.org/hprc-v2.1-mc-grch38.GRCh38.walks.gz',
)

const refName = 'GRCh38#0#chr6'
const { graph, kept, fragments, departures } = await reader.cut(
  refName,
  { start: 160_614_798, end: 160_647_758 },
  { keep: walksKept(['HG00097', 'HG00133'], refName) },
)
const tables = graph.tables(kept, fragments, departures)
```

- `WalkReader.open(url)` reads the walk file at `url` and its directory at `url`
  plus `.idx`, put before any query string or fragment, by HTTP range requests
  through generic-filehandle2's `RemoteFile`. Its second argument opens each
  location another way: `url => new LocalFile(url)` for a local file
- `WalkReader.fromFilehandles({ data, directory })` takes the two handles
  already open: `data` the `.walks.gz`, `directory` its `.idx`. Any object with
  a generic-filehandle2 handle's `read`, `readFile`, `stat` and `source` will do
- A file gfa-to-tabix 4.0 or later wrote; a format 6 `.walks.bed.gz` has no
  directory, and the reader says to rebuild it. jbrowse.org hosts the HPRC
  release 2.1 graph's walk file in format 7,
  `https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38.GRCh38.walks.gz`
- `tables` lists the nodes with their lengths and reference positions, the links
  between them, each walk's steps, and every departure from the reference over
  the window ([Departures](#departures)). `graphTablesGFA` in
  [@jbrowse/bandage-core](https://github.com/GMOD/bandage-core) writes the
  nodes, links and walks as GFA
- A node the reference passes more than once, as in a collapsed repeat, takes
  the position of its pass over the window, and each walk's step on it the pass
  its row names
- `reader.header()` lists the reference sample and every haplotype with rows,
  and the tile size, `tile`, read from the walk file's header unit
- `reader.referenceSequences()` lists the sequences a region can be on, in full
  PanSN: `GRCh38#0#chr6`

## Changes in 6

graph-walks 6 reads walk format 7 only (gfa-to-tabix's ADR 0010): a series of
gzip units with a unit table per sequence and a directory, not a Tabix file.

- Open a walk file with `WalkReader.open(url)` or
  `WalkReader.fromFilehandles({ data, directory })`, in place of
  `new WalkReader(new TabixIndexedFile(...))`. `cut`, `cutRegion`, `header`,
  `referenceSequences` and `sequenceFor` keep their calls, and `cut` its return
  shape
- A cut reads the tiles and departure blocks overlapping the window and its
  context, by extent, and nothing before them: no chunk before the window, no
  `maxnode` lookback, no `#long` widening, no reach on the departure read, no
  `~departure` or `~cs` sequences. An attach fetch reads the tiles a record's
  `ownTiles` names, for the record's own steps alone. Units whose bytes touch
  come in one range request
- Since a cut reads less, a fragment's `trimmedBefore` and `trimmedAfter`, the
  bp of its walk the cut read and kept out, can be smaller than graph-walks 5
  gave for the same window
- An abridged fragment's `abridged.detail` names its detail unit: the sequence,
  the tile's nominal start and `name`, `tile:detail`
- The byte budgets count the table's exact unit lengths, so they hold on every
  source
- New: `reader.units(refName)`, a sequence's unit table;
  `reader.sequenceLength(refName)`, where its reference ends;
  `reader.bytesFor(refName, region)`, the bytes a cut there reads;
  `reader.referenceNodes(refName, region)`, every pass of the reference over it;
  `WalkFile`, `WalkUnit`, `WalkFilehandle` and `WALK_FORMAT`
- Removed: `queryStart`, `chunkQueryStart`, `lookbackChunks`, `walkQueryStart`,
  `LongPass`, and the header's `chunk`, `maxNode` and `long`; `byteBudgetError`
  takes a synchronous `bytesTo`. `TabixSource` now types only the segment index
- Links are written the way round with more `+` ends, as a walk along the
  reference takes them, and with the smaller first handle on a tie: a walk
  `>7>3` gives `L 7 + 3 +`, where graph-walks 5 gave `L 3 - 7 -`
- Fixed: two departure records whose leave pieces lie in adjacent tiles both get
  their own steps; graph-walks 5.0.2 dropped one's

## Changes in 5

graph-walks 5 reads walk format 6 only and changes the API:

- `cut` returns `departures`, every departure from the reference over the window
  once, in place of `skips`, `open` and `across`; a fragment runs from its
  walk's first step on the window to its last, with no follow rounds, crossings
  or stop kinds. `WalkGraph.tables` takes the departures
- Removed: `CROSSING_BP`, `WalkChunk`, `WalkCrossing`, `WalkGap`, `WalkSkip`,
  `WalkStop`, `WalkStopTable`; `WalkFragment`, `WalkRow`, `GraphTables` and
  `HaplotypeDepartures` lost fields, and `keptStepsError` takes the departures
  as its second argument
- Budgets default to `WALK_BUDGETS`, so a cut that read without limit now throws
  `RegionTooLargeError` or `DepartureReadError` past them; pass `Infinity` to
  lift one
- A departure states how many times the reference passes each end's node, and
  `placementChoice` says where its position was the writer's choice
- Structural forms state nothing a member's own walks contradict: a rare state
  folds only into one it agrees with, a no call never folds, a member no common
  form fits goes to `other` (key `*`), and a change at a site wider than 1 kb
  keys on its place. `formContradictions` checks a result against the members'
  departures

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
| `stepBudget`          | The most steps a cut keeps, departures' own steps included; 4,000,000 by default    |
| `byteBudget`          | The most compressed bytes the walk and attach reads fetch; 8 MB by default          |
| `departureByteBudget` | The most compressed bytes the departure read fetches; 2 MB by default               |
| `attachBp`            | The most bp of own sequence a departure gets its steps attached for; default 65,536 |
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
- The cut checks the node limit as it decodes each tile, against the nodes it
  keeps whatever else it reads, and stops decoding at the first tile over it
- A window past `departureByteBudget` throws a `DepartureReadError`, which is
  not a `RegionTooLargeError`: a departure's record is filed over its whole
  stretch, so zooming in does not lower that read. A collapsed repeat that
  haplotypes hop between the copies of is the usual cause
- `WALK_BUDGETS` holds the default budgets and `attachBp`, what a browser window
  holds; a caller that wants no limit passes `Infinity`
- The byte budgets bound what a cut fetches: the unit table gives the exact
  compressed bytes of the window's tiles and of its departure blocks before
  either is fetched, and the cut refuses either past its budget. The tiles
  fetched for records' own steps count against `byteBudget` with the window's,
  and a cut they would take past it throws the same `RegionTooLargeError` before
  fetching any of them, so a window attaches every eligible record's own steps
  or none. The detail records `expand` names count only against `stepBudget`,
  once read; the writer moves paths out of any tile whose rows would compress
  past its header's `abridgebytes`
- A narrow window reads whole tiles, so a cut reads several times the steps it
  keeps. The rows a cut decodes from the window's tiles, and from the detail
  records `expand` names once they are decoded, may hold up to 8 times
  `stepBudget` steps (`READ_STEPS_PER_KEPT`), a step along the reference
  counting a quarter; past that the cut throws a `RegionTooLargeError` before it
  cuts. The rows of tiles fetched for records' own steps are not counted there,
  and the own steps they give count against `stepBudget` with the rest

## Fragments

Each of `fragments` is one walk's stretch through the cut: `name`, `fragStart`
(the W line's start), `hapStart`/`hapEnd`, and its steps, from the walk's first
step on the drawn stretch of the reference, the window and its context, to its
last. `trimmedBefore` and `trimmedAfter` count the bp of the walk's run past
either end of the fragment, so they depend on the tiles the cut read. They
include own steps the cut attaches to a departure there, so they say how far the
run goes, not what lies there: [Gaps](#gaps) does.

- Fragments come from the tiles overlapping the drawn stretch alone. A tile the
  cut fetches for a record's own steps (see [Departures](#departures)) gives
  that record its steps and nothing else, so no fragment holds a step from it
- A walk's pieces in those tiles join where their piece indices are consecutive.
  Where a walk goes into a tile the cut did not read and comes back, it has a
  fragment either side, each ending at its last step on the stretch, and its
  departure records say where it leaves and rejoins
- A walk's fragments therefore depend on the window, its context and the file
  alone, not on which other walks the cut keeps
- `graph.tables(kept, fragments, departures)` keeps `trimmedBefore` and
  `trimmedAfter` as columns of `walks`

## Gaps

`gaps` states what lies on each walk's contig between its fragments, past its
first and last, and in its departures outside them, as
[ADR 0001](docs/adr/0001-a-cut-states-what-lies-between-fragments.md) decides.
Each has `name`, `edge` (0 between two fragments, 1 a row's start, 2 its end, 3
a departure outside the row), `start`/`end` on the contig, the `fragment` it
follows (or precedes, for a row's start), and `segments` in contig order whose
`bp` sum to `end - start`:

- `steps`: own steps the cut attached to a departure, with their `ids` and
  `rev`; these nodes are in the cut
- `read`: steps of the walk's runs over the window's tiles that the cut kept
  out, such as reference steps off the drawn stretch or own steps of a departure
  over `attachBp`
- `abridged`: an abridged run of those tiles that is not a fragment
- `own`: a record's own sequence the cut did not read, with its `missingBp`, in
  no stated order, and `abridged` set where the record's leave or rejoin piece
  is an abridged run
- `unread`: pieces the cut did not read, with the reference breakpoints the cut
  knows `before` and `after` them

A row end runs over the bp of the run past the fragment, and over the whole of
any record that reaches it. Tiles fetched for a record's own steps give `steps`
alone, and the cut fetches them for every record within `attachBp` or refuses
the window, so which records get `steps` follows from each record alone.
`graph.tables(kept, fragments, departures, gaps)` writes them as a `gaps` table,
each stretch's walk by its fragment's index in `walks`, and each step as a node
index.

## Departures

A departure is where a walk's consecutive steps on the reference are not
adjacent there: it leaves the reference at one breakpoint and rejoins it at
another, with `ownBp` of its own sequence between. `cut` returns each departure
once:

- Every departure record whose line overlaps the window and its context, from
  the departure blocks there. gfa-to-tabix writes one for every departure
  between two pieces of a walk. A record comes back even where no reference node
  lies under the window, as in a `#hole`
- Every departure inside one piece of a walk, which the cut computes from the
  piece's steps, whose stretch, from the base beside its leave to the base
  beside its rejoin, overlaps the drawn stretch: the reference's nodes over the
  window and its context, whole, widened by every pass of a walk's step over
  them. A walk's jump inside one piece far from the window, even inside its
  fragment, is not returned

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
- `leavePasses` and `rejoinPasses` say how many times the reference passes each
  end's node. `placementChoice(d)` is true where either is above 1: the writer
  chose which pass the step takes, so the departure is its choice, not a jump
  the walk makes, as where a walk follows one copy of a collapsed repeat the
  reference passes several times. A node's `passes` in `graph.nodes`, and the
  tables' `nodes.passes`, say the same of each step on the reference, from the
  ref rows the cut read; a node the cut places off the reference has 0, though
  it may lie on the reference outside those rows; the header's `revisited`
  counts the nodes passed more than once, 0 where every departure in the file is
  a fact of the graph
- A departure of at most `attachBp` of own sequence carries its own steps in
  `steps`: the walk's steps between its leave and its rejoin, from the end of
  the leave piece and the start of the rejoin piece. A step of 0 bp is an own
  step, so a record with `ownBp` 0, or no more than its `missingBp`, can carry
  one; own sequence between the contig's W lines has no steps
- They come from the pieces the cut read. For a record whose leave or rejoin
  piece lies outside the read, the cut fetches the tiles its `ownTiles` column
  names that it has not read for that walk, and takes from them that record's
  own steps alone. It fetches them for every such record, or throws a
  `RegionTooLargeError` where they would take the cut past `byteBudget`, so
  whether a record gets its steps does not depend on the other walks the cut
  keeps, as gfa-to-tabix ADR 0009 asks
- Own steps in an abridged piece come only for a walk `expand` names, which
  reads its detail unit
- `fromSteps` is true for a departure computed from a piece's steps
- Those steps' nodes are kept, so the cut draws a short insertion or replacement
  at the window's edge whole, and counts them against `stepBudget`
- `graph.tables(kept, fragments, departures)` keeps them as the `departures`
  table: `names`, `hapStarts`, `ownBp`, `missingBp`, `leavePos`,
  `leaveReversed`, `leavePasses`, `rejoinPos`, `rejoinReversed`, `rejoinPasses`,
  `fromSteps`, and each one's own steps as node indices at `offsets[k]` to
  `offsets[k + 1]` of `steps`, with `reversed`

## Abridged walks

- In a tile over its budget, gfa-to-tabix moves the costliest paths' rows to
  detail records and leaves one `abridged` row per moved piece. A satellite
  array in a pggb graph is the usual case
- The cut makes a separate fragment for each such row, with an `abridged` field:
  the piece's `steps` and `bp`, and the `detail` unit that holds it with the
  path's other moved pieces in the tile: its sequence `refName`, the tile's
  nominal `start`, and its `name`, `tile:detail`. The fragment's two steps are
  the piece's first and last handles on the reference, `hapStart`/`hapEnd` are
  exactly where the stretch between them lies on its contig, and `steps` counts
  all the piece's steps, its own before and after that stretch included
- `graph.tables(kept, fragments, departures)` keeps them as `abridgedSteps` and
  `abridgedBp` columns, NaN for every other walk
- With `expand` naming a walk, the cut reads its detail records and cuts it in
  full, as a file without abridging would
- `walkDepartures` nets an unexpanded piece over its whole stretch, `bp` against
  the reference between its handles, which may reach past the window; an
  expanded one is trimmed to the drawn stretch like any walk, so the two reads
  differ over a piece straddling the drawn stretch's edge

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
- `formsCut` takes `byteBudget`, and a region whose tiles the table puts past it
  throws `RegionTooLargeError` before any read. The budget covers the region's
  tiles, not the detail records an expanded cut reads, which only the step
  budget bounds
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

| Option          | Effect                                                                                                                    |
| --------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `window`        | The region; marks must touch it, and `auto` reads its width                                                               |
| `floorBp`       | The smallest change a form tells apart; `auto` gives 50 bp, 1 kb or 10 kb by width                                        |
| `minCarriers`   | The fewest haplotypes a state or form needs; 1% of them, at least one, by default                                         |
| `joinBp`        | Departures of one haplotype nearer than this on the reference join; default 50                                            |
| `exclude`       | Samples to leave out, such as a second reference carried as a haplotype; none by default                                  |
| `foldRareForms` | Fold states and forms under `minCarriers` into the nearest common one they agree with, else `~` or `other`; on by default |

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

- By default a cut keeps each walk from its first step on the drawn stretch, the
  window and its context, to its last, whatever it visits between, and follows
  nothing past that stretch, so one haplotype's fragment never depends on which
  others are cut
- `contained: true` splits a walk wherever it steps onto the reference off the
  drawn stretch. Each piece between those steps that steps on the stretch is
  kept from its first such step to its last, so a detour that comes back stays
  whole, with `context` bp of the walk past either end
- A fragment counts the bp of its walk the cut read past either end as
  `trimmedBefore` and `trimmedAfter`
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
  file's sequence names

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

- A tile holds the reference's path through it once, each way a haplotype leaves
  that path once, and every other path's piece as the list of ways it takes;
  gfa-to-tabix's
  [walks-format.md](https://github.com/GMOD/gfa-to-tabix/blob/main/docs/walks-format.md)
  describes the rows, and its ADR 0010 the file
- At open the reader reads the directory, refusing any format but 7, then the
  header unit and the sequence's unit table at once, and refuses a directory
  whose first line differs from the header's `#walks` line or whose `bytes`
  differs from the data's length. Over HTTP the length comes from Content-Range,
  so a server must expose that header to cross-origin readers for the check to
  hold; where the reader cannot learn the length it skips the check. Every call
  that reads the file, `referenceSequences`, `units` and `referenceNodes` among
  them, waits on this check
- It refuses, with a `WalkFormatError`, a unit or table that does not hold what
  lists it: one that is not whole gzip members, one whose text ends partway
  through a line, a unit whose `#unit` line is not its table entry's, and a row
  of another sequence or of a kind its unit does not hold. A unit of a kind the
  format does not give it never reads
- A cut selects the `tile` and `departures` units whose `start end` overlap the
  window and its context, skipping any of length 0, merges the byte ranges that
  touch, and fetches each as one single-range request. It gunzips each unit,
  decodes each tile once, and keeps the departure lines overlapping the window
  and its context, each line once
- A departure inside one piece comes from the piece's steps, and one between
  pieces from its record, so each comes back once. A departure's own steps are
  attached where its own bp are under `attachBp`, from its leave piece's steps
  after the leave and its rejoin piece's before the rejoin. Where the read lacks
  either piece, the cut fetches the tiles the record's `ownTiles` names, every
  one wanted in one more read counted against `byteBudget`, and takes from them
  the records' own steps alone; a tile fetched again decodes only the walks it
  was not decoded for. Expanding an abridged piece reads its `detail` unit

## License

MIT
