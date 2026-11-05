# Structural forms from walks

- **In:** a `gfa-to-tabix --walks` file and a reference window
- **Out:** the window's forms: groups of haplotypes that make the same
  structural changes, each change with its kind, exact REF span, ALT length and
  carriers
- **Where:** `walkDepartures.ts` and `structuralForms.ts`. The graphgenomeviewer
  plugin (GPL) and jb2hubs call them
- **Measured on:** the published `hprc-v2.1-mc-grch38-f5.GRCh38.walks.bed.gz`,
  Node 24 on a laptop, 462 haplotypes once GRCh38 and CHM13 are left out

## Sites come from the walks

Each haplotype's departures from the reference path R give the sites. The
alternatives fall short:

| Source                               | Verdict                                                                                                                                                                       |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| bubbles file                         | Top-level bubbles only, with no per-haplotype allele; an rGFA's comes from `gfatools` on the SV graph, whose segment ids differ from the base-level graph the walk file holds |
| MinigraphBubbleAdapter / GetSubgraph | One feature per bubble with no genotypes                                                                                                                                      |
| **walks against R**                  | **Chosen**                                                                                                                                                                    |

`walkDepartures` reads each haplotype off the cut's departures and its
fragments' ends, in chains:

- A chain is a fragment with the departures leaving its last step, rejoining its
  first or lying inside it, fragments that abut on the contig, or a departure
  touching no fragment, such as a contig that jumps over the window or rejoins
  far outside it, which the window's chunks never held. Each chain reads in its
  own orientation: one that reads R mostly reversed reads the other way round (5
  of 3,740 fragments measured)
- Between a departure's rejoin and the next leave the path follows R. A rejoin
  that is not where the last stretch along R ended closes a departure. Its REF
  span runs from that end to the rejoin, and the path's own bp between, sequence
  the graph lacks between a contig's W lines included, give its ALT length
- A stretch of R walked against it inside a forward chain counts as ALT and as
  inverted bp (AMY1: 317 haplotypes)
- A step back along R is a tandem copy. The departure is the zero-width
  junction, and the bp walked again count as ALT, as `repeatedBp`, with their
  span as `copied`
- R that two chains of one haplotype reach, such as a contig's second visit to
  the window or a second contig over the same locus, is no call there
- An abridged fragment's bp lie between the R its handles reach, so the
  departure holding it is exact in net, and only placement inside it is lost.
  Each abridged row is one piece, from its hapOffset to hapOffset + bp, so a
  piece of another chunk between two moved pieces is a fragment of its own and
  the departures between them are records. Before gfa-to-tabix wrote a row per
  piece, one row stood for a path's moved pieces in a chunk with the bp of those
  pieces alone, and a reader had to add what lay between them; leaving that out
  put AMY1's unexpanded HG03710#1 29 kb short: +167,764 bp against +196,894
  expanded
- A departure leaving or rejoining inside a piece keeps all its own bp, and the
  piece's net leaves out those of its steps: its trailing steps after the leave,
  its leading steps before the rejoin. The piece's handles on R lie outside
  every departure's own sequence, so that is all they share. A piece runs along
  R or against it as both its handles do, so one against R in a chain read
  forward is inverted ALT, as its steps would be
- The reader this replaced kept, of a departure with an end inside a piece, only
  the bp the graph lacks and the own steps the cut held, which it attaches only
  up to `ATTACH_BP`, and joined two pieces a departure left and rejoined inside
  into one net. Own steps in no piece of the cut, such as a W line holding none
  of R or a step on another reference sequence, counted nowhere, and a piece's
  trailing own steps went into its net even where the departure rejoined
  reversed and the walks read in full close nothing. On MC chr22 that put 11
  haplotypes in 29 windows off their expanded net: HG03453#1 at
  chr22:18,219,008-18,284,544 read +227,316 bp unexpanded against +2,382
  expanded, HG00329#2 256,134 against 268,533. Unexpanded and expanded now give
  all 463 haplotypes the same net and coverage in each of those 29 windows, the
  28 of the forms arbitration and that one

Walks give exact lengths. Each haplotype's summed departures across the truth
table's span match its assembly in 2,779 of 2,780 scored cases, median gap 0 bp;
the HPRC sidecar's equivalent is 99% with a median of 36 bp. The one miss,
NA19909#1 at AMY1 (walk −94,529 bp, assembly −101 bp), came from a second contig
over the AMY copies whose departures added to the first's; R that two walks of
one haplotype cover is now no call.

## Grouping rules

| Rule            | Decision                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `stateKey`      | Size change to 2 significant figures, for grouping only: one LPA copy number differs by a few bp (22,201 vs 22,204). A shift keys on its size, and on its net change where that reaches the floor and the size's second figure, 1 kb for a 94 kb shift: under that the net is within the size's own rounding                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| Detection floor | 50 bp. At a site no wider than 1 kb it applies to each haplotype's net change across the site, its plain pieces summed, since walks split a repeat's change into pieces: at rand29's VNTR the walk net matches the assembly for 408 of 408 haplotypes, the pieces that pass 50 bp alone for 235. A wider site keeps its changes apart: summed, the haplotypes without UGT2B17's 117 kb deletion added up small indels into 7 forms of −50 to −56 bp                                                                                                                                                                                                                                                                                                                                                                                                    |
| `minCarriers`   | max(1, round(n / 100)): 5 of 462                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `changeKind`    | On exact REF/ALT with no rounding slack: an insertion where REF is under 50 bp, a deletion where ALT is, otherwise a replacement                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Join            | Departures of one haplotype less than 50 bp of R apart join, keeping their net. The join replaces the sidecar's split-event merging: at GSTM1 walks give 306 carriers of an exact −18,445 bp deletion and 156 with no change, where the sidecar had 6 forms, two of them "shift" forms holding 82 haplotypes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `pairShifts`    | On exact marks: a deletion of 1 kb or more and an insertion or duplication within 10% of its size and no further than twice it. At AMY1, 191 haplotypes pair a −94,204 bp deletion with a +94,142 bp insertion 18 kb away; 164 of them are among the sidecar's 171 "another route" haplotypes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Nesting         | Marks that `minCarriers` haplotypes share, two marks shared where each covers at least half of the other on R, are common, and overlapping common marks form one site, whose span is theirs; a haplotype's state there is only its own marks, so a parent deletion and a nested child never both count. UGT2B17's forms held 230/229/3 haplotypes from walks and 232/229/1 from the sidecar under the earlier union of every mark, not remeasured. A rarer mark states itself at each site it overlaps without widening it: at chr22:20.3-20.4 Mb one HG01175#1 departure over 108.7 kb of R had merged two sites of 461 and 460 haplotypes into one of 137 states. A summed change reaching past its site does not widen it either: at chr22:50.78 Mb a widened site held 1 kb insertions of 3 haplotypes whose state was read over the narrower span |
| Place           | At a site wider than `SUMMED_SITE_BP` each change keys on its place as well as its size, `+14000@42125912`: the first start among the marks the half-cover relation joins it to. Spans that cover half of each other share a place; insertions share one only at one point. CYP2D6's 13.7 kb insertion has place 42,125,912 in 38 haplotypes and 42,134,095 in 4 (MC chr22 42.10-42.15 Mb, 1 kb floor); keyed on size alone they were one state, drawn at the first place for all 42. A summed site keys on its net change alone, so a repeat's change keeps one state wherever the walks put it                                                                                                                                                                                                                                                       |
| Rare states     | A state under `minCarriers` folds into the common state nearest in net change among those it agrees with: change by change at a site keyed on place, net and largest change at a summed site, each pair under `agrees` (both under the floor, or apart by under the floor, a tenth of the larger, or the last figure the key gives the change's size). A change no common state agrees with is `~`, rarer. R's structure and a no call stay themselves. `foldRareForms: false` keeps every state, as a scorer wants every member's own sizes                                                                                                                                                                                                                                                                                                           |
| Rare forms      | A haplotype whose form has under `minCarriers` joins the largest common form whose state at every site is its own or one its changes agree with, or `other`, the last form, key `*`, which states nothing of its members. Without it LPA at 50 bp has 48 single-haplotype rows                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Edge sites      | A site that most haplotypes in the window have no call at is `edge` and left out of the key. A haplotype's call reads its coverage over the site's part inside the window, so a departure followed past the edge does not make the others no-calls                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| No calls        | A haplotype whose walks leave a gap in a site's part inside the window has no call there, even with a mark at the site, as a walk running off R at AMY1's insertion junction has. R that two of its walks cover, such as a second contig, and a junction between two walks, such as a record break whose bases the graph lacks, are gaps too                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| Shared sites    | A site where every form but `other` is in one state, a call, tells none apart: it is out of the key and `formTag` leaves its marks out                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |

### What a form states

A form states, at every site that is not an edge, one state of its members: its
marks there, a no call, or a rarer change. `formContradictions` checks that
statement against each member's own departures, independently of how the sites
and states were built:

- `net`: the member's net change at the site, its joined departures touching it
  at a summed site and its marks touching it at a wider one, agrees with the
  state's marks
- `size`: so does its largest change
- `place`: at a wider site, its marks touching each of the state's marks agree
  with the state's marks touching it
- `call`: the form states a call only where the member's walks span the site's
  part inside the window, and `hidden` a no call only where they do not
- `rarer`: a rarer change only where the member makes one
- `unsited`: every mark of the member lies in a site

`other` and edge sites state nothing, so nothing there can contradict.
structuralForms meets it by construction: a state folds, and a haplotype joins a
form, only through the same `agrees`, and what agrees with nothing is `~` or
`other`.

Members a result contradicts, summed over windows, edf4acd's nearest fold
against these rules:

| Windows                                 | Floor, fold     | Nearest fold | These rules |
| --------------------------------------- | --------------- | ------------ | ----------- |
| 28 MC chr22 windows, format 6           | auto, folded    | 315          | 0           |
| 28 MC chr22 windows, format 6           | 50 bp, folded   | 1,849        | 0           |
| 28 MC chr22 windows, format 6           | 50 bp, unfolded | 663          | 0           |
| jb2hubs gate, 19 loci and 40 random, f6 | auto, folded    | 417          | 0           |
| jb2hubs gate, 19 loci and 40 random, f6 | 50 bp, folded   | 1,073        | 0           |
| jb2hubs gate, 19 loci and 40 random, f6 | 50 bp, unfolded | 206          | 0           |

The cost is members that no common form fits: at the auto floor 70 in `other`
and 371 in a form with a rarer change, of 25,887 placed over the gate's 59
windows; at 50 bp folded, 561 and 502. They claim no size there. The gate's
claims, 50 bp unfolded, scored by `scoreForms.py --complete-claims`: 6,981 of
7,003 at the loci and 13,013 of 13,013 in the random windows, against 6,931 of
6,952 and 12,946 of 12,946 with the nearest fold, which scored 51 and 67 rare
carriers apart and missed 4 and 18 of them. The 28 windows are the `W` lines of
gfa-to-tabix's `mc22.forms-compare-review.txt`.

The nearest fold this replaced folded any rare state to the common state nearest
in net, a rare form to the form nearest in net, and a rare no call to the
majority. That kept the rows few and named, and at LPA, AMY1 and MHC raised the
prototype's claims over a majority fold (LPA 381 to 395 of 406, MHC 367 to 407);
the part of that gain that is honest stays, since a state within tolerance still
folds, and the rest becomes `~` or `other`, which claims no size.

### Mark kinds

| Kind                                   | Meaning                                                                    |
| -------------------------------------- | -------------------------------------------------------------------------- |
| `insertion`, `deletion`, `replacement` | `changeKind` on the exact REF and ALT                                      |
| `inversion`                            | Net change under the floor, inverted bp at or over it                      |
| `duplication`                          | A step back along R; `repeatedBp` and `copied` give the copy               |
| `shift`                                | A paired deletion and insertion or duplication; `shiftedBp` gives the size |
| `net`                                  | An abridged walk left unexpanded: its net change somewhere in its span     |

A replacement whose smaller side is under 5% of its larger one
(`LOPSIDED_RATIO`) keeps kind `replacement`, its REF span, ALT length and delta,
so a scorer sees the exact change. It gains `lopsided: 'deletion'` or
`'insertion'`, which a caller can colour as a loss or gain, and `markTag` writes
it as one: `del 6.46kb->94bp` for C4's 6,460 bp to 94 bp.

Duplication is a kind of its own because core JBrowse already has the term
(`DUP:TANDEM` folds to `DUP`, and VcfFeature maps it to `tandem_duplication`),
so the plugin can colour them as a gain. AMY1 has 231 back-steps of 50 bp or
more over 160 haplotypes, 19-65 kb each (median 48.6 kb); LPA and MHC have none,
since their extra copies show as insertions. Whether a copy reads as a back-step
depends on how the graph collapsed the repeat, so the key treats a duplication
like an insertion of its size. Moving its REF onto the copied span instead would
join neighbouring sites: AMY1 goes from 18 to 23 forms (ARI 0.915).

## Labels

`markTag` writes a mark as a short token, `formTag` a form as its largest token
with ` +N` for the rest, and `describeForm` every token, largest first; the
README's table lists them. `sizeTag` writes 3 significant figures, so
`del 6.46kb->94bp` gives both sides of C4's change. The earlier net label read
"−6.4 kb" and dropped the 94 bp left over.

Row labels at the `auto` floor, 462 haplotypes once CHM13 is left out:

| Window | Haplotypes | `formTag`                  | Note                                                      |
| ------ | ---------- | -------------------------- | --------------------------------------------------------- |
| C4     | 217        | `del 6.46kb->94bp`         | Lopsided: 94 bp left of 6,460                             |
| C4     | 95         | `shift 32.7kb`             | `onlyShifts`                                              |
| C4     | 20         | `repl 6.46kb->26.5kb`      | Over 5% on both sides                                     |
| CYP2D6 | 42         | `ins 178bp->13.9kb`        | Lopsided: 178 bp replaced                                 |
| MHC    | 26         | `del 87.5kb->920bp +2`     | Lopsided                                                  |
| MHC    | 31         | `repl 8.95kb->50.9kb +10`  | `describeForm` lists all 11 tokens                        |
| AMY1   | 175        | `shift 94.2kb`             | H3r, `onlyShifts`                                         |
| AMY1   | 16         | `shift 94.2kb net +2.81kb` | `onlyShifts`; the net passes 1 kb, the size's last figure |
| AMY1   | 61         | `dup 94.1kb`               |                                                           |

`onlyShifts` marks a form whose changes outside the shared sites are all shifts,
with a call at every site: R's length, with a stretch aligned elsewhere. AMY1's
two such forms hold 191 of 460 placed haplotypes, so a display can set them
beside the reference row. A form with a no call or a rarer state is not one,
since its walks leave part of the window unread.

## Display floor

`floorBp: 'auto'` takes the window's width / 100, about 10 px, snapped to 50 bp
under 25 kb, 1 kb from 25 to 300 kb and 10 kb above. Forms that become identical
at a higher floor merge.

Forms, every one with 5 or more haplotypes once rare forms fold, measured under
the nearest fold and not remeasured:

| Window      | 50 bp | 500 bp | 1 kb | 5 kb |
| ----------- | ----- | ------ | ---- | ---- |
| GSTM1 35 kb | 2     | 2      | 2    | 2    |
| C4 70 kb    | 11    | 10     | 10   | 10   |
| LPA 130 kb  | 26    | 30     | 30   | 30   |
| AMY1 150 kb | 13    | 14     | 14   | 14   |
| MHC 90 kb   | 29    | 21     | 14   | 9    |

At gene scale a fixed 1 kb floor leaves the track blank: in HP
72,040,000-72,045,000, a 300 bp deletion in 5 haplotypes falls to one "no
change" row, and a rand1 tile with a 51-55 bp insertion does the same. The auto
rule gives that 5 kb window 50 bp and the five windows above 1 kb.

## Abridged walks

`formsCut` expands abridged walks, fewest steps first, under `FORMS_STEP_BUDGET`
(16 M steps), and leaves any past it as `net` marks.

- Net-only is a poor default: at MHC 304 haplotypes get `net` marks spanning a
  whole 60-65 kb chunk, and at LPA only 153 haplotypes can be scored at all (406
  expanded)
- An expanded MHC cut keeps 9.4 M steps, past `WALK_BUDGETS`' 4 M, so `formsCut`
  has its own budget, which its cuts keep to in place of the reader's default
- Expansion costs a second cut over the same chunks, since only the first lists
  the abridged walks and their steps
- `partial: false` reads every abridged walk or none. Where reading them all
  would take the cut past `stepBudget`, `formsCut` throws
  `ExpansionTooLargeError` after the plain cut, before the second. Partial
  expansion is a poor middle: MHC at a 4 M budget costs 3.5 s and still leaves
  191 haplotypes net-only, worse than the sidecar
- A gfa-to-tabix `dp` tag listing departures on abridged rows would remove the
  expansion. Readers skip unknown tags (`walkChunk.ts`), so it needs no format
  bump, but the HPRC files would need a rebuild: 39 min and 13 GB on 24 cores,
  plus the 64 GB base-level GFA

`expand` false against true, measured by jb2hubs, wall time:

| Window | Plain cut | All expanded | RSS               |
| ------ | --------- | ------------ | ----------------- |
| MHC    | 1.7-2.4 s | 9.5 s        | 0.41 → 1.5-1.7 GB |
| KIR    | 3.6 s     | 7.3 s        |                   |
| LPA    | 1.9 s     | 4.8 s        |                   |
| AMY1   | 1.7 s     | 3.1 s        |                   |

Expanded cuts before parallel detail reads, and `formsCut` with 16 at once
(graph-walks after 4.0.0), cold / warm:

| Window      | Plain cut     | Abridged walks | Expanded, serial    | `formsCut`    | Steps expanded |
| ----------- | ------------- | -------------- | ------------------- | ------------- | -------------- |
| GSTM1 35 kb | 0.85 / 0.24 s | 0              | (none)              | 0.46 / 0.31 s | 0.4 M          |
| LPA 130 kb  | 0.82 / 0.74 s | 291            | 8.4-11.6 / 2.1 s    | 4.8 / 3.5 s   | 6.1 M          |
| AMY1 150 kb | 0.72 / 0.60 s | 132            | 3.8-4.3 / 1.2-1.7 s | 2.5 / 2.1 s   | 4.0 M          |
| MHC 90 kb   | 1.0 / 1.0 s   | 519            | 19-25 / 9-15 s      | 8.2 / 7.4 s   | 9.4 M          |

`structuralForms` takes under 0.1 s in every window, so the cut dominates.
`formsCut` spends about 1 s over a single expanded cut on its first, plain cut.

### What a `dp` tag must list

structuralForms reads more than the 50 bp departures, so a `dp` tag of those
alone would change the forms:

- Departures already joined at `joinBp` (50 bp), as `join` joins them, so each
  entry is one change with its exact net
- The nonzero pieces under 50 bp inside each change of 50 bp or more. At a site
  up to `SUMMED_SITE_BP` (1 kb) wide, structuralForms sums a haplotype's plain
  pieces, which lifts rand29's VNTR from 235 haplotypes matching the assembly to
  408 of 408
- Gate it on forms from `dp` equalling the expanded forms at MHC, KIR, LPA, AMY1
  and the jb2hubs gate windows

## Validation

Through `formsCut` and `structuralForms` against the published file, expanded,
scored against the assemblies (`truth.tsv`) and the sidecar (`forms-v3.jsonl`):

| Window | Walk net = assembly | Forms 50 bp / 1 kb | ARI vs sidecar, 50 bp | Sidecar forms |
| ------ | ------------------- | ------------------ | --------------------- | ------------- |
| GSTM1  | 405/405             | 2 / 2              | 0.86                  | 6             |
| C4     | 407/407             | 11 / 10            | 0.95                  | 15            |
| LPA    | 406/406             | 26 / 30            | 0.37                  | 179           |
| MHC    | 407/407             | 29 / 14            | 0.67                  | 120           |

MHC's net is scored on a cut widened to the truth table's span,
chr6:32,476,308-32,610,000. The ARI is low where the sidecar splits LPA into 179
forms and MHC into 120, most of them single haplotypes.

The jb2hubs gate (`printWalkForms.ts`, `scoreForms.py --complete-claims`) over
19 loci and 40 random 50 kb windows, at 50 bp with `foldRareForms: false`, on
the published format 5 file before the nearest fold was replaced (format 6
figures under [What a form states](#what-a-form-states)):

- Claims within 10% or 500 bp of the assembly: 6,910 of 6,931 at the loci and
  12,946 of 12,946 in the random windows. The sidecar sizes 5,528 and 12,860
- Per site, 50 bp to 1 kb, members both score: walks 78.8% against the sidecar's
  69.9% at the loci, 96.1% against 92.7% in the random windows
- Of the 21 misses, 20 are at KIR, where a walk and the assembly's flanks put a
  16 kb or 31 kb stretch on opposite sides of a span's edge among the KIR
  paralogs. The other, NA18508#1 at FCGR, walks one contig 82 kb longer than the
  flanks measure

## Consumers

- The plugin's `WalkFormsAdapter` feeds core `LinearMultiSampleVariantDisplay`,
  one row per form, pinned to a window, with the floor in its menu. It caches
  `HaplotypeDepartures` by walk file and pinned window, since a new floor
  regroups them without a new cut
- jb2hubs moves the HPRC and cattle pages to walk forms; cattle needs a walk
  file built from its P lines. Arabidopsis keeps the sidecar until its graph has
  paths
