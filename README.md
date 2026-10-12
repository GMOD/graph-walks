# @jbrowse/graph-walks

Cuts a window out of a pangenome graph stored as
[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) walk files, by range
requests.

- **In:** a `gfa-to-tabix --walks` file set and a region on its reference
- **Out:** the reference over the region, every haplotype's walk through it, and
  the nodes and links they visit, as typed arrays
- Reads format 1 (`.walks.bed.gz` with `.nodes.bed.gz` and `.links.bed.gz`) and
  format 2 (`.walks.bed.gz` alone), each with its `.tbi`
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
import {
  WalkReader,
  panSNHaplotype,
  walkNameFilter,
} from '@jbrowse/graph-walks'
import { RemoteFile } from 'generic-filehandle2'

const file = (url: string) =>
  new TabixIndexedFile({
    filehandle: new RemoteFile(url),
    tbiFilehandle: new RemoteFile(`${url}.tbi`),
  })
const prefix = 'https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38.GRCh38'
const reader = new WalkReader({
  walks: file(`${prefix}.walks.bed.gz`),
  nodes: file(`${prefix}.nodes.bed.gz`),
})

const refName = 'GRCh38#0#chr6'
const { graph, kept, fragments } = await reader.cut(
  refName,
  { start: 160_614_798, end: 160_647_758 },
  { keep: walkNameFilter(['HG00097', 'HG00133'], panSNHaplotype(refName)) },
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
| `links`      | `false` leaves the link file unread                                       |
| `sequences`  | Keep each node's bases in `graph.sequences`, from `--sequences` files     |
| `signal`     | An `AbortSignal` that stops the reads                                     |

A cut over a budget throws a `RegionTooLargeError`. Its `fitsBp` holds about how
wide a window would fit, and JBrowse shows it as a zoom-in notice.

## How a cut works

- A cut reads every chunk under the region, and the chunks before it as far as
  the header's `maxnode` reaches, since a node is filed under the chunk holding
  its start
- A haplotype that leaves the reference near the region is followed to where it
  rejoins, and the reference is cut out as far as those rejoins reach
- A walk the read cuts off goes on in the chunk its row names (`pv:`/`nx:`
  tags); the cut reads on to that chunk, in up to three rounds

## License

MIT
