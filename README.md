# @jbrowse/graph-walks

Cuts a window out of a pangenome graph stored as
[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) walk files, by range
requests.

- **In:** a `gfa-to-tabix --walks` file, `<prefix>.<sample>.walks.bed.gz` with
  its `.tbi`, and a region on its reference
- **Out:** the reference over the region, every haplotype's walk through it, and
  the nodes and links they visit, as typed arrays
- Reads walk format 3, which gfa-to-tabix writes from 1.1. The reader refuses
  any other format with a `WalkFormatError` asking for a rebuild
- Has no runtime dependencies: it reads through any object with
  [@gmod/tabix](https://github.com/GMOD/tabix-js)'s `TabixIndexedFile` methods

The JBrowse graph plugin, BandageJS and the sequenceTubeMap fork read walk files
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

## Cut options

| Option       | Effect                                                                    |
| ------------ | ------------------------------------------------------------------------- |
| `keep`       | Which walks to decode, by `sample#haplotype#contig`; every one by default |
| `context`    | bp of reference kept past each end of the region; default 1000            |
| `stepBudget` | The most steps a cut keeps; no limit by default                           |
| `byteBudget` | The most compressed bytes a cut fetches; no limit by default              |
| `links`      | `false` lists no links, for a caller that draws no edges                  |
| `sequences`  | Keep each node's bases in `graph.sequences`, from a `--sequences` file    |
| `signal`     | An `AbortSignal` that stops the reads                                     |

A cut over a budget throws a `RegionTooLargeError`. Its `fitsBp` holds about how
wide a window would fit, and JBrowse shows it as a zoom-in notice.
`WALK_BUDGETS` holds the step and byte budgets the JBrowse plugin and BandageJS
cut by.

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

## Naming the window

- `reader.sequenceFor(sample, refName)` gives the walk file's sequence for a
  contig on a sample: `chr6` on `GRCh38` is `GRCh38#0#chr6`, and a pggb file's
  `chm13#chr1` is `chr1` on `chm13`
- `walksKept(haplotypes, refName)` is the `keep` filter for a cut on `refName`:
  its own walk and those of the haplotypes named
- `refNameLookup` and `resolveRefName` resolve a contig the same way against any
  tabix file's sequence names

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
