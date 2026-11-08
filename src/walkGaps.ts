import type { WalkDeparture } from './walkDeparture.ts'
import type { WalkFragment, WalkRun } from './walkRows.ts'

export type GapKind = 'steps' | 'read' | 'abridged' | 'own' | 'unread'

export const GAP_KINDS: GapKind[] = [
  'steps',
  'read',
  'abridged',
  'own',
  'unread',
]

/**
 * One segment of a walk's stretch, by ADR 0001: `steps` the cut attached as
 * own sequence, with their ids; `read`, steps the cut read and kept out;
 * `abridged`, an abridged run that is not a fragment; `own`, a record's
 * unread middle, with its `missingBp` and whether an end piece was abridged;
 * `unread`, pieces the cut did not read, between the breakpoints on the
 * reference it knows `before` and `after` them
 */
export interface GapSegment {
  kind: GapKind
  bp: number
  missingBp: number
  abridged: boolean
  before: number
  after: number
  ids: number[]
  rev: number[]
}

/**
 * A stretch of a walk's contig and its segments: `edge` 0 between two
 * fragments, 1 a row's start, 2 a row's end, 3 a departure outside the row;
 * `fragment` the one it follows, or precedes for a row's start
 */
export interface WalkGap {
  name: string
  edge: 0 | 1 | 2 | 3
  start: number
  end: number
  fragment?: WalkFragment
  segments: GapSegment[]
}

interface Span {
  start: number
  end: number
}

interface Step extends Span {
  id: number
  rev: number
  attached: boolean
}

interface Anchor extends Span {
  run: WalkRun
  i: number
}

interface Placed extends Span {
  d: WalkDeparture
}

export interface GapInputs {
  runs: WalkRun[]
  attach: WalkRun[]
  fragments: WalkFragment[]
  departures: WalkDeparture[]
  refName: string
  before: (run: WalkRun) => Float64Array
  length: (id: number) => number
  attached: Map<WalkRun, Set<number>>
  spans: { run: WalkRun; first: number; last: number }[]
  breakpoint: (run: WalkRun, i: number, leave: boolean) => { pos: number }
}

function groupByName<T extends { name: string }>(items: T[], refName: string) {
  const out = new Map<string, T[]>()
  for (const item of items) {
    if (item.name !== refName) {
      const list = out.get(item.name)
      if (list) {
        list.push(item)
      } else {
        out.set(item.name, [item])
      }
    }
  }
  return out
}

// the first index in sorted `spans` whose end passes `at`, or touches it
function firstEnding<T extends Span>(spans: T[], at: number, touching = false) {
  let lo = 0
  let hi = spans.length
  while (lo < hi) {
    const mid = (lo + hi) >> 1
    const end = spans[mid]!.end
    if (touching ? end >= at : end > at) {
      hi = mid
    } else {
      lo = mid + 1
    }
  }
  return lo
}

// the first index in `spans`, sorted by start, that starts at or past `at`
function firstStarting<T extends Span>(spans: T[], at: number) {
  let lo = 0
  let hi = spans.length
  while (lo < hi) {
    const mid = (lo + hi) >> 1
    if (spans[mid]!.start >= at) {
      hi = mid
    } else {
      lo = mid + 1
    }
  }
  return lo
}

const runEnd = (run: WalkRun, bp: Float64Array) =>
  run.hapOffset + (run.abridged?.bp ?? bp[run.ids.length]!)

/**
 * Each walk's stretches between and past its fragments, and its departures
 * outside them, split into segments that partition each stretch, as
 * docs/adr/0001-a-cut-states-what-lies-between-fragments.md decides
 */
export function walkGaps({
  runs,
  attach,
  fragments,
  departures,
  refName,
  before,
  length,
  attached,
  breakpoint,
  spans,
}: GapInputs) {
  const inFragment = new Map<WalkRun, [number, number][]>()
  for (const { run, first, last } of spans) {
    inFragment.set(run, [...(inFragment.get(run) ?? []), [first, last]])
  }
  const fragmentStep = (run: WalkRun, i: number) =>
    inFragment.get(run)?.some(([first, last]) => i >= first && i <= last) ??
    false
  const fragmentAbridged = new Set(
    fragments.flatMap(f => (f.abridged ? [f.abridged] : [])),
  )
  const runsBy = groupByName(runs, refName)
  const attachBy = groupByName(attach, refName)
  const fragmentsBy = groupByName(fragments, refName)
  const departuresBy = groupByName(departures, refName)
  const names = [
    ...new Set([...fragmentsBy.keys(), ...departuresBy.keys()]),
  ].sort()

  const gaps: WalkGap[] = []
  for (const name of names) {
    const walkRuns = [...(runsBy.get(name) ?? [])].sort(
      (a, b) => a.hapOffset - b.hapOffset,
    )
    const stepRuns = walkRuns.filter(run => !run.abridged)
    const abridgedSpans = walkRuns.flatMap(run =>
      run.abridged
        ? [
            {
              start: run.hapOffset,
              end: run.hapOffset + run.abridged.bp,
              fragment: fragmentAbridged.has(run.abridged),
            },
          ]
        : [],
    )
    const attachSteps: Step[] = []
    for (const run of attachBy.get(name) ?? []) {
      const marks = attached.get(run)
      if (marks) {
        const bp = before(run)
        for (const i of marks) {
          attachSteps.push({
            start: run.hapOffset + bp[i]!,
            end: run.hapOffset + bp[i + 1]!,
            id: run.ids[i]!,
            rev: run.rev[i]!,
            attached: true,
          })
        }
      }
    }
    const placed: Placed[] = (departuresBy.get(name) ?? [])
      .map(d => ({ d, start: d.hapStart, end: d.hapStart + d.ownBp }))
      .sort((a, b) => a.start - b.start || a.end - b.end)
    const walkFragments = [...(fragmentsBy.get(name) ?? [])].sort(
      (a, b) => a.hapStart - b.hapStart,
    )

    // the reference steps of the walk's runs, and an abridged run's two
    // handles, built once a stretch needs a breakpoint
    let anchors: Anchor[] | undefined
    const anchorsOf = () => {
      anchors ??= walkRuns
        .flatMap(run => {
          const bp = before(run)
          return run.abridged
            ? [
                {
                  run,
                  i: 0,
                  start: run.hapOffset,
                  end: run.hapOffset + length(run.ids[0]!),
                },
                {
                  run,
                  i: 1,
                  start: runEnd(run, bp) - length(run.ids[1]!),
                  end: runEnd(run, bp),
                },
              ]
            : [...run.ids.keys()]
                .filter(i => !Number.isNaN(run.refStarts[i]))
                .map(i => ({
                  run,
                  i,
                  start: run.hapOffset + bp[i]!,
                  end: run.hapOffset + bp[i + 1]!,
                }))
        })
        .sort((a, b) => a.start - b.start)
      return anchors
    }
    const breakpointsAround = (start: number, end: number) => {
      const all = anchorsOf()
      const back = firstEnding(all, start) - 1
      const lastAnchor = back >= 0 ? all[back] : undefined
      const nextAnchor = all[firstStarting(all, end)]
      const backPlaced = firstEnding(placed, start) - 1
      const lastPlaced = backPlaced >= 0 ? placed[backPlaced] : undefined
      const nextPlaced = placed[firstStarting(placed, end)]
      const beforeFromPlaced =
        lastPlaced && (!lastAnchor || lastPlaced.end > lastAnchor.end)
      const afterFromPlaced =
        nextPlaced && (!nextAnchor || nextPlaced.start < nextAnchor.start)
      return {
        before: beforeFromPlaced
          ? lastPlaced.d.rejoinPos
          : lastAnchor
            ? breakpoint(lastAnchor.run, lastAnchor.i, true).pos
            : Number.NaN,
        after: afterFromPlaced
          ? nextPlaced.d.leavePos
          : nextAnchor
            ? breakpoint(nextAnchor.run, nextAnchor.i, false).pos
            : Number.NaN,
      }
    }

    // a stretch's read steps, a 0 bp one included where it lies on either
    // edge
    const stepsIn = (start: number, end: number) => {
      const out: Step[] = []
      for (const run of stepRuns) {
        const bp = before(run)
        if (run.hapOffset <= end && runEnd(run, bp) >= start) {
          const marks = attached.get(run)
          let lo = 0
          let hi = run.ids.length
          while (lo < hi) {
            const mid = (lo + hi) >> 1
            if (run.hapOffset + bp[mid + 1]! >= start) {
              hi = mid
            } else {
              lo = mid + 1
            }
          }
          for (
            let i = lo;
            i < run.ids.length && run.hapOffset + bp[i]! <= end;
            i++
          ) {
            const s = run.hapOffset + bp[i]!
            const e = run.hapOffset + bp[i + 1]!
            const isAttached = marks?.has(i) ?? false
            if (
              !fragmentStep(run, i) &&
              ((e > start && s < end) || (s === e && isAttached))
            ) {
              out.push({
                start: s,
                end: e,
                id: run.ids[i]!,
                rev: run.rev[i]!,
                attached: isAttached,
              })
            }
          }
        }
      }
      for (const x of attachSteps) {
        if (
          (x.end > start && x.start < end) ||
          (x.start === x.end && x.start >= start && x.start <= end)
        ) {
          out.push(x)
        }
      }
      return out.sort((a, b) => a.start - b.start || a.end - b.end)
    }

    // how far a row's end reaches over the departures that touch it
    const reach = (start: number, end: number, forward: boolean) => {
      let at = forward ? end : start
      for (const x of forward ? placed : [...placed].reverse()) {
        if (forward && x.end > start && x.start <= at && x.end > at) {
          at = x.end
        } else if (!forward && x.start < end && x.end >= at && x.start < at) {
          at = x.start
        }
      }
      return at
    }

    const stretches: Omit<WalkGap, 'segments'>[] = []
    const first = walkFragments[0]
    const last = walkFragments.at(-1)
    const rowStart = first
      ? reach(first.hapStart - first.trimmedBefore, first.hapStart, false)
      : Infinity
    const rowEnd = last
      ? reach(last.hapEnd, last.hapEnd + last.trimmedAfter, true)
      : -Infinity
    if (first && rowStart < first.hapStart) {
      stretches.push({
        name,
        edge: 1,
        start: rowStart,
        end: first.hapStart,
        fragment: first,
      })
    }
    walkFragments.forEach((f, k) => {
      const next = walkFragments[k + 1]
      if (next) {
        stretches.push({
          name,
          edge: 0,
          start: f.hapEnd,
          end: next.hapStart,
          fragment: f,
        })
      }
    })
    if (last && rowEnd > last.hapEnd) {
      stretches.push({
        name,
        edge: 2,
        start: last.hapEnd,
        end: rowEnd,
        fragment: last,
      })
    }
    for (const x of placed) {
      const outside =
        x.end > x.start && (!first || x.end <= rowStart || x.start >= rowEnd)
      if (outside) {
        stretches.push({ name, edge: 3, start: x.start, end: x.end })
      }
    }
    stretches.sort((a, b) => a.start - b.start || a.edge - b.edge)

    for (const stretch of stretches) {
      const { start, end } = stretch
      const inSteps = stepsIn(start, end)
      const sized = inSteps.filter(x => x.end > x.start)
      const zeroSteps = inSteps.filter(x => x.end === x.start && x.attached)
      const inAbridged = abridgedSpans.filter(
        x => !x.fragment && x.end > start && x.start < end,
      )
      const from = firstEnding(placed, start, true)
      const inPlaced: Placed[] = []
      for (let k = from; k < placed.length && placed[k]!.start <= end; k++) {
        inPlaced.push(placed[k]!)
      }
      const sizedPlaced = inPlaced.filter(x => x.end > start && x.start < end)
      const zeroPlaced = inPlaced.filter(
        x => x.end === x.start && x.start > start && x.start < end,
      )
      const points = new Set([start, end])
      for (const x of [
        ...sized,
        ...inAbridged,
        ...sizedPlaced,
        ...zeroPlaced,
      ]) {
        for (const p of [x.start, x.end]) {
          if (p > start && p < end) {
            points.add(p)
          }
        }
      }
      const sorted = [...points].sort((a, b) => a - b)
      const segments: GapSegment[] = []
      let lastOwner: WalkDeparture | undefined
      const add = (
        kind: GapKind,
        bp: number,
        more: Partial<GapSegment> = {},
        owner?: WalkDeparture,
      ) => {
        const prev = segments.at(-1)
        if (
          bp > 0 &&
          prev?.kind === kind &&
          prev.bp > 0 &&
          (kind !== 'own' || lastOwner === owner)
        ) {
          prev.bp += bp
          prev.after = more.after ?? prev.after
          for (const [j, id] of (more.ids ?? []).entries()) {
            prev.ids.push(id)
            prev.rev.push(more.rev![j]!)
          }
        } else {
          segments.push({
            kind,
            bp,
            missingBp: 0,
            abridged: false,
            before: Number.NaN,
            after: Number.NaN,
            ids: [],
            rev: [],
            ...more,
          })
          lastOwner = owner
        }
      }
      let p = 0
      let z = 0
      const flushZeros = (at: number) => {
        for (; z < zeroSteps.length && zeroSteps[z]!.start <= at; z++) {
          const x = zeroSteps[z]!
          add('steps', 0, { ids: [x.id], rev: [x.rev] })
        }
        for (const x of zeroPlaced) {
          if (x.start === at) {
            add('own', 0, { missingBp: x.d.missingBp }, x.d)
          }
        }
      }
      for (let k = 0; k + 1 < sorted.length; k++) {
        const a = sorted[k]!
        const b = sorted[k + 1]!
        flushZeros(a)
        while (p < sized.length && sized[p]!.end <= a) {
          p++
        }
        const candidate = sized[p]
        const step =
          candidate && candidate.start <= a && candidate.end >= b
            ? candidate
            : undefined
        const piece = inAbridged.find(x => x.start <= a && x.end >= b)
        const own = sizedPlaced.find(x => x.start <= a && x.end >= b)
        const bp = b - a
        if (step) {
          add(
            step.attached ? 'steps' : 'read',
            bp,
            step.attached ? { ids: [step.id], rev: [step.rev] } : {},
          )
        } else if (piece) {
          add('abridged', bp)
        } else if (own) {
          add(
            'own',
            bp,
            {
              missingBp: own.d.missingBp,
              abridged: abridgedSpans.some(
                x => x.end === own.start || x.start === own.end,
              ),
            },
            own.d,
          )
        } else {
          add('unread', bp, breakpointsAround(a, b))
        }
      }
      flushZeros(end)
      gaps.push({ ...stretch, segments })
    }
  }
  return gaps
}
