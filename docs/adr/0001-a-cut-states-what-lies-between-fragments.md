# 0001: A cut states what lies between a walk's fragments

## Status

Proposed, 2026-11-08

## Context

A walk row in @jbrowse/bandage-core draws a walk's fragments end to end along
its contig, and the stretches between them as gaps. Until now a consumer worked
out what a gap holds from numbers the cut returns for other jobs:

- `trimmedBefore` and `trimmedAfter` count the bp of a run that lie past its
  fragment. That can be own steps the cut attaches to a departure, own steps of
  a departure over `attachBp` (read, never attached), or reference steps off the
  drawn stretch.
- A record's `steps` joins the steps after its leave and the steps before its
  rejoin, with no boundary between them.
- Abridged runs report trimmed counts of 0 and never take attached steps.
- Computed departures attach their steps as one run.

At the amylase locus, HG00097#1's 29,895 bp of own steps, attached to a record
and kept in the cut, were drawn as "outside the extracted subgraph", because the
consumer read `trimmedAfter` as sequence the cut did not keep. Eight reviews of
rules that rebuild the layout from these counts each found another case they get
wrong: several records in one gap, abridged neighbours, partial attaches,
records over `attachBp`, contained cuts, and row ends. `walkCut` already holds
what those rules try to rebuild.

## Decision

The cut returns, for each walk, stretches of its contig, each split into
segments in contig order. Every segment is a fact the cut holds, and a stretch's
segments partition it exactly. A consumer draws them as given.

A walk's stretches are:

- **Between two fragments:** from one fragment's `hapEnd` to the next one's
  `hapStart`. A stretch of 0 bp, where two pieces abut, is still a row, holding
  only any 0 bp steps that lie there.
- **A row's start and end:** from a fragment's edge over the bp of its run the
  cut read past it, `trimmedBefore` or `trimmedAfter`. The stretch is extended
  to cover the whole own stretch of any record it reaches. A row end of 0 bp is
  left out.
- **A record outside the row:** the own stretch of a record that lies wholly
  outside the walk's other stretches, and every record of a walk with no
  fragment. The graph draws these records' attached nodes, so the cut states
  them too. Each consumer decides how to show them.

Each base of a stretch falls in one segment, by the first rule that applies:

- **`steps`**: the base lies on a step the cut attached as own sequence, either
  to a record within `attachBp` or to a computed departure the cut returns
  within `attachBp`. The segment carries those of the departure's steps that lie
  in the stretch, as node indices and orientations; these nodes are in the cut.
  A 0 bp step makes a 0 bp segment at its position, at either edge of the
  stretch too.
- **`read`**: the base lies on any other step of the walk's runs over the
  window's tiles. Examples are a reference step off the drawn stretch, or an own
  step of a departure over `attachBp`. The cut read it and kept it out. Runs the
  cut fetched only for a record's own steps (C132) contribute `steps` alone,
  never `read` or `abridged`, so those never depend on which records the cut
  fetched for.
- **`abridged`**: the base lies in the stretch of an abridged run that is not a
  fragment. The cut did not load its steps; expanding the walk does.
- **`own`**: the base lies in a record's own stretch,
  `[hapStart, hapStart + ownBp)`. The segment is the record's whole unread
  middle. It carries the record's `missingBp`, the part of its own sequence
  between the contig's W lines that the graph lacks. It also carries whether the
  record's leave or rejoin piece was an abridged run among the runs over the
  window's tiles. The file does not state where within that middle the missing
  bp lie, or the order of its steps on other reference sequences or on W lines
  with no step on the reference (walk format C80, C87, ADR 0009 of
  gfa-to-tabix), so the segment states the amounts and no order.
- **`unread`**: any other base. The walk lies in pieces the cut did not read.
  Those pieces may hold departures the cut did not return, because it returns
  only those over the drawn stretch (C135), so the segment states no more than
  its bp and the breakpoints the cut knows either side:
  - `before`: the leave breakpoint of the last step on the reference in the
    walk's runs before it, or the rejoin of the last record before it, whichever
    is nearer.
  - `after`: the rejoin breakpoint of the next such step, or the leave of the
    next record.

  Breakpoints are positions on the reference, as `breakpoint()` takes them, so a
  reversed step's is its node's start. They say where the walk was last and next
  seen, not where the unread pieces lie.

A record of 0 bp inside a stretch makes a 0 bp `own` segment at its position, so
the walk is placed there; one at a stretch's edge, and one outside the row,
makes nothing. Segments of the same kind and more than 0 bp that abut are
merged, except `own` segments of different records.

## Shape

`GraphTables` gains a `gaps` table, one row per stretch, in walk then contig
order:

- `names`: the walk's name, so a walk with no fragment has stretches too
- `walks`: the walk-table index of the fragment the stretch follows, or of the
  fragment it precedes for a row's start; -1 for a record outside the row
- `edges`: 0 for a stretch between two fragments, 1 for a row's start, 2 for a
  row's end, 3 for a record outside the row
- `starts`, `ends`: the stretch on the contig
- `segmentOffsets`: segment rows `segmentOffsets[g]` up to
  `segmentOffsets[g + 1]`

And a segment table:

- `kinds`: 0 `steps`, 1 `read`, 2 `abridged`, 3 `own`, 4 `unread`
- `bp`
- `missingBp`, `abridged`: for `own`, 0 otherwise
- `before`, `after`: for `unread`, NaN otherwise or where the cut knows no
  breakpoint on that side
- `stepOffsets`, `steps`, `reversed`: for `steps`, its node indices and
  orientations

`WalkFragment`'s `trimmedBefore` and `trimmedAfter` stay, documented as the bp
of the run read past the fragment. They no longer describe what those bp are.

A consumer that rebuilds walks as paths and filters them, as bandage-core's
Graph path does, carries each walk's stretches on its path, not by walk-table
index.

## Consequences

- A consumer draws a gap without inferring from counts. In bandage-core,
  `splitGap`'s trimmed arithmetic and its per-case rules go.
- HG00097#1's row draws its 131 attached steps, so it shares the copy HG00099#1
  walks, then an `own` segment of 19,156 bp carrying 19,156 missing bp.
- Every `read` and `unread` segment depends only on the walk's runs over the
  window's tiles and on the records returned. `steps` and `own` also depend on
  `attachBp`, and on one more thing: the attach-fetch cap (`attachFetches`,
  first come) decides which records' far pieces the cut fetches. A record whose
  piece the cap leaves unfetched shows its unread bp as `own`, which is true,
  but the same event can look different for two haplotypes. That conflicts with
  gfa-to-tabix ADR 0009's "never from budget order" and is filed separately.
- A contained cut, whose fragments end `context` bp into a run, follows the same
  rules. Steps between two fragments of one run are read, so they fall in
  `steps` or `read`.

## Alternatives considered

- **Per-departure columns** (`leaveSteps` and a status for each side), with the
  consumer laying out gaps. This keeps the layout as consumer inference over
  trimmed counts, which is what failed review after review.
- **Keeping a W line's off-reference end steps in its fragment.** It cannot tell
  a line's end from a jump into a tile the cut did not read without consulting
  records, and it changes what a fragment is for every consumer.
- **Relabelling in the consumer only.** It draws nothing false, but a row then
  disagrees with the graph, which does draw the kept nodes.
