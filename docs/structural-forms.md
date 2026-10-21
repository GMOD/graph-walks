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

`walkDepartures` steps along R in each walk:

- A fragment that reads R mostly reversed reads the other way round (5 of 3,740
  fragments measured)
- A forward R step that does not abut the last one closes a departure. Its REF
  span runs from the last R step's end to this one's start, and the steps
  between give its ALT length
- A reversed R step inside a forward walk counts as ALT and as inverted bp
  (AMY1: 317 haplotypes)
- A step back along R is a tandem copy. The departure is the zero-width
  junction, and the bp walked again count as ALT, as `repeatedBp`, with their
  span as `copied`
- Off-R steps before the first R step or after the last leave the walk open: the
  haplotype has no call there
- Fragments of one path that abut on the haplotype read as one walk. An abridged
  fragment's bp then lie between the R steps either side of it, so the departure
  holding it is exact in net, and only placement inside it is lost
- An abridged row's pieces need not be consecutive. Where the path steps onto
  another chunk's R and back, the pieces between lie in that chunk, and their bp
  join the abridged departure's. Leaving them out put AMY1's unexpanded
  HG03710#1 29 kb short: +167,764 bp against +196,894 expanded
- An abridged fragment at either end of a walk has no R step beyond it, so its
  departure stops at its anchors, which is exact only where the walk enters or
  leaves R there. At AMY1, 129 of 132 unexpanded abridged haplotypes match their
  expanded net and coverage exactly

Walks give exact lengths. Each haplotype's summed departures across the truth
table's span match its assembly in 2,779 of 2,780 scored cases, median gap 0 bp;
the HPRC sidecar's equivalent is 99% with a median of 36 bp. The one miss,
NA19909#1 at AMY1 (walk −94,529 bp, assembly −101 bp), has no explanation yet.

## Grouping rules

| Rule            | Decision                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `stateKey`      | Size change to 2 significant figures, for grouping only: one LPA copy number differs by a few bp (22,201 vs 22,204). A shift keys on its size, and on its net change where that reaches the floor and the size's second figure, 1 kb for a 94 kb shift: under that the net is within the size's own rounding                                                                                                                                        |
| Detection floor | 50 bp. At a site no wider than 1 kb it applies to each haplotype's net change across the site, its plain pieces summed, since walks split a repeat's change into pieces: at rand29's VNTR the walk net matches the assembly for 408 of 408 haplotypes, the pieces that pass 50 bp alone for 235. A wider site keeps its changes apart: summed, the haplotypes without UGT2B17's 117 kb deletion added up small indels into 7 forms of −50 to −56 bp |
| `minCarriers`   | max(1, round(n / 100)): 5 of 462                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `changeKind`    | On exact REF/ALT with no rounding slack: an insertion where REF is under 50 bp, a deletion where ALT is, otherwise a replacement                                                                                                                                                                                                                                                                                                                    |
| Join            | Departures of one haplotype less than 50 bp of R apart join, keeping their net. The join replaces the sidecar's split-event merging: at GSTM1 walks give 306 carriers of an exact −18,445 bp deletion and 156 with no change, where the sidecar had 6 forms, two of them "shift" forms holding 82 haplotypes                                                                                                                                        |
| `pairShifts`    | On exact marks: a deletion of 1 kb or more and an insertion or duplication within 10% of its size and no further than twice it. At AMY1, 191 haplotypes pair a −94,204 bp deletion with a +94,142 bp insertion 18 kb away; 164 of them are among the sidecar's 171 "another route" haplotypes                                                                                                                                                       |
| Nesting         | Overlapping marks form one site, and a haplotype's state there is only its own marks, so a parent deletion and a nested child never both count. UGT2B17's forms hold 230/229/3 haplotypes from walks and 232/229/1 from the sidecar                                                                                                                                                                                                                 |
| Rare states     | Fold to the common state nearest in net change. The sidecar's restored parents are dropped                                                                                                                                                                                                                                                                                                                                                          |
| Rare forms      | A form under `minCarriers` folds into the form nearest in net. Without it LPA at 50 bp has 48 single-haplotype rows                                                                                                                                                                                                                                                                                                                                 |
| Edge sites      | A site that most haplotypes in the window have no call at is `edge` and left out of the key. A haplotype's call reads its coverage over the site's part inside the window, so a departure followed past the edge does not make the others no-calls                                                                                                                                                                                                  |
| No calls        | A haplotype whose walks leave a gap in a site's part inside the window has no call there, even with a mark at the site, as an open walk at AMY1's insertion junction has                                                                                                                                                                                                                                                                            |
| Shared sites    | A site where no second state reaches `minCarriers` is not in the key, and every form carries its majority state's marks                                                                                                                                                                                                                                                                                                                             |

### Folding to the nearest state

The prototype, before rare forms folded, at a floor of 50 bp; AMY1 and MHC on
cuts widened to the truth table's span:

| Window | Claims, majority fold | Claims, nearest fold | Folded haplotypes, majority → nearest |
| ------ | --------------------- | -------------------- | ------------------------------------- |
| LPA    | 381/406               | 395/406              | 35 → 49 of 60                         |
| AMY1   | 239/274               | 245/274              | 4 → 10 of 23                          |
| MHC    | 367/407               | 407/407              | 3 → 43 of 43                          |

Folding rare forms as well costs claims at 50 bp, since the form nearest in net
can differ by one KIV-2 copy: LPA 395 → 385, MHC 407 → 406. At 1 kb, the floor
`auto` gives both windows, LPA keeps 395.

### Mark kinds

| Kind                                   | Meaning                                                                    |
| -------------------------------------- | -------------------------------------------------------------------------- |
| `insertion`, `deletion`, `replacement` | `changeKind` on the exact REF and ALT                                      |
| `inversion`                            | Net change under the floor, inverted bp at or over it                      |
| `duplication`                          | A step back along R; `repeatedBp` and `copied` give the copy               |
| `shift`                                | A paired deletion and insertion or duplication; `shiftedBp` gives the size |
| `net`                                  | An abridged walk left unexpanded: its net change somewhere in its span     |

Duplication is a kind of its own because core JBrowse already has the term
(`DUP:TANDEM` folds to `DUP`, and VcfFeature maps it to `tandem_duplication`),
so the plugin can colour them as a gain. AMY1 has 231 back-steps of 50 bp or
more over 160 haplotypes, 19-65 kb each (median 48.6 kb); LPA and MHC have none,
since their extra copies show as insertions. Whether a copy reads as a back-step
depends on how the graph collapsed the repeat, so the key treats a duplication
like an insertion of its size. Moving its REF onto the copied span instead would
join neighbouring sites: AMY1 goes from 18 to 23 forms (ARI 0.915).

## Display floor

`floorBp: 'auto'` takes the window's width / 100, about 10 px, snapped to 50 bp
under 25 kb, 1 kb from 25 to 300 kb and 10 kb above. Forms that become identical
at a higher floor merge.

Forms, every one with 5 or more haplotypes once rare forms fold:

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
  has its own budget
- Expansion costs a second cut over the same chunks, since only the first lists
  the abridged walks and their steps
- A gfa-to-tabix `dp` tag listing departures of 50 bp or more on abridged rows
  would remove the expansion. Readers skip unknown tags (`walkChunk.ts`), so it
  needs no format bump, but the HPRC files would need a rebuild: 39 min and 13
  GB on 24 cores, plus the 64 GB base-level GFA

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

## Consumers

- The plugin's `WalkFormsAdapter` feeds core `LinearMultiSampleVariantDisplay`,
  one row per form, pinned to a window, with the floor in its menu. It caches
  `HaplotypeDepartures` by walk file and pinned window, since a new floor
  regroups them without a new cut
- jb2hubs moves the HPRC and cattle pages to walk forms; cattle needs a walk
  file built from its P lines. Arabidopsis keeps the sidecar until its graph has
  paths
