import { panSNMatchesPrefix } from './util.ts'

import type { Departure, HaplotypeDepartures } from './walkDepartures.ts'

/**
 * What a haplotype's change does to R's bases over its REF span: inserts at a
 * point, deletes, replaces them, inverts them, or walks a stretch of R again
 * (`duplication`). A `shift` pairs a deletion with an insertion of about its
 * size near it, and `net` is an abridged walk's change somewhere in its span.
 */
export type MarkKind =
  | 'insertion'
  | 'deletion'
  | 'replacement'
  | 'inversion'
  | 'shift'
  | 'duplication'
  | 'net'

// `delta` is altLength minus the REF span. A duplication states the bp it
// walks again and their span on R, and a shift the size of its deletion.
export interface Mark {
  kind: MarkKind
  refStart: number
  refEnd: number
  altLength: number
  delta: number
  repeatedBp?: number
  copied?: [number, number]
  shiftedBp?: number
}

// A mark in a structuralForms result: `site` indexes its `sites`, and `state`
// is the key of the haplotype's state there
export interface SiteMark extends Mark {
  site: number
  state: string
}

/**
 * One state of a site: its key, the median net change of its carriers, and
 * the marks of the carrier with that change. `carriers` are the haplotypes
 * whose own marks give the state, before any folding. `0` is R's structure.
 */
export interface SiteState {
  key: string
  delta: number
  marks: SiteMark[]
  carriers: string[]
}

/**
 * Overlapping marks make one site. `covered` counts the placed haplotypes
 * with a call there: walks spanning the site's part inside the window, so a
 * mark beside a gap in them is no call. An `edge` site is one most haplotypes in the window have no call
 * at, which no form key holds. A `shared` site tells no forms apart, since no
 * second state there reaches minCarriers, and every form carries its majority
 * state's marks. Every other site is informative.
 */
export interface Site {
  refStart: number
  refEnd: number
  edge: boolean
  shared: boolean
  covered: number
  states: SiteState[]
}

/**
 * Haplotypes with one state at every informative site. `states` holds a
 * SiteState key per informative site, `.` (NO_CALL_STATE) where the form's
 * walks do not span it, and `~` (RARER_STATE) for a rare change where most
 * haplotypes have no call; `key` joins them. `marks` are its states' marks in
 * R order, with those of the state every form shares at a site that tells
 * none apart, and `net` sums them.
 */
export interface Form {
  members: string[]
  key: string
  states: string[]
  marks: SiteMark[]
  net: number
}

export interface FormsOptions {
  // the smallest change a form tells apart; `auto`, the default, follows the
  // window's width (autoFloor)
  floorBp?: number | 'auto'
  // departures of one haplotype nearer than this on R join
  joinBp?: number
  minCarriers?: number
  window?: { start: number; end: number }
  // samples or `sample#haplotype`s to leave out, such as a second reference
  // carried as a haplotype
  exclude?: string[]
  // fold forms under minCarriers into the nearest common form; true by default
  foldRareForms?: boolean
}

/**
 * The forms a window's haplotypes fall into. `informative` indexes the sites
 * where a second state has minCarriers. `placed` counts the haplotypes in a
 * form, and `unplaced` lists those with no call in the window and no form.
 * `rareCarriers` carry a state or a form too rare to stand alone, folded into
 * the nearest common one, and `abridged` have a `net` mark.
 */
export interface StructuralForms {
  sites: Site[]
  informative: number[]
  forms: Form[]
  placed: number
  unplaced: string[]
  rareCarriers: string[]
  abridged: string[]
  minCarriers: number
  floorBp: number
}

// The smallest change detected, and the side of a change under which it
// inserts or deletes rather than replaces
export const STRUCTURAL_BP = 50

export const FLOOR_PRESETS = [50, 1000, 10_000] as const

const SHIFT_MIN_BP = 1000
const SHIFT_TOLERANCE = 0.1

// A form's state at a site where it keeps R's structure, where its walks have
// no call, and where it makes a change too rare to stand alone while most
// haplotypes have no call
export const REFERENCE_STATE = '0'
export const NO_CALL_STATE = '.'
export const RARER_STATE = '~'

/**
 * The floor for a window `widthBp` wide: about 10 px at 1,000 px, snapped to a
 * preset, so 50 bp under 25 kb, 1 kb to 300 kb and 10 kb past it
 */
export function autoFloor(widthBp: number) {
  const px = widthBp / 100
  const [fine, medium, coarse] = FLOOR_PRESETS
  return px < 250 ? fine : px < 3000 ? medium : coarse
}

export const minCarriersFor = (haplotypes: number) =>
  Math.max(1, Math.round(haplotypes / 100))

// Under STRUCTURAL_BP on one side a change inserts or deletes; otherwise it
// replaces R's bases
export function changeKind(refLength: number, altLength: number) {
  return altLength > refLength && refLength < STRUCTURAL_BP
    ? 'insertion'
    : altLength < refLength && altLength < STRUCTURAL_BP
      ? 'deletion'
      : 'replacement'
}

const twoFigures = (n: number) => {
  const unit = 10 ** (Math.floor(Math.log10(n)) - 1)
  return Math.round(Math.round(n / unit) * unit)
}

/**
 * A mark's state, its size change to two significant figures: one repeat
 * unit is one state wherever the walk's boundaries fall (22,201 and 22,204 bp
 * are both `+22000`), and two units another. An inversion keys on its span,
 * a shift on its size, and a net change apart from a placed one.
 */
export function stateKey(m: Mark) {
  const sign = m.delta > 0 ? '+' : '-'
  switch (m.kind) {
    case 'inversion':
      return `v${twoFigures(Math.max(Math.abs(m.delta), m.refEnd - m.refStart, 1))}`
    case 'shift':
      return `s${twoFigures(m.shiftedBp!)}`
    case 'net':
      return `n${sign}${twoFigures(Math.abs(m.delta))}`
    default:
      return `${sign}${twoFigures(Math.abs(m.delta))}`
  }
}

interface Joined {
  refStart: number
  refEnd: number
  delta: number
  invertedBp: number
  repeatedBp: number
  copied?: [number, number]
  abridged?: true
}

const deltaOf = (d: Departure) => d.altLength - (d.refEnd - d.refStart)

// Departures nearer than `joinBp` on R joined, keeping their net change
function join(departures: Departure[], joinBp: number) {
  const out: Joined[] = []
  for (const d of departures) {
    const last = out.at(-1)
    if (
      last &&
      !last.abridged &&
      !d.abridged &&
      d.refStart - last.refEnd < joinBp
    ) {
      last.refEnd = Math.max(last.refEnd, d.refEnd)
      last.delta += deltaOf(d)
      last.invertedBp += d.invertedBp
      last.repeatedBp += d.repeatedBp
      if (d.copied) {
        last.copied = last.copied
          ? [
              Math.min(last.copied[0], d.copied[0]),
              Math.max(last.copied[1], d.copied[1]),
            ]
          : d.copied
      }
    } else {
      out.push({
        refStart: d.refStart,
        refEnd: d.refEnd,
        delta: deltaOf(d),
        invertedBp: d.invertedBp,
        repeatedBp: d.repeatedBp,
        ...(d.copied ? { copied: d.copied } : {}),
        ...(d.abridged ? { abridged: true as const } : {}),
      })
    }
  }
  return out
}

function markOf(d: Joined, floorBp: number): Mark | undefined {
  const refLength = d.refEnd - d.refStart
  const altLength = Math.max(0, refLength + d.delta)
  const size = Math.abs(d.delta)
  const kind: MarkKind | undefined = d.abridged
    ? size >= floorBp
      ? 'net'
      : undefined
    : size < floorBp
      ? d.invertedBp >= floorBp
        ? 'inversion'
        : undefined
      : d.repeatedBp >= STRUCTURAL_BP && d.delta > 0
        ? 'duplication'
        : changeKind(refLength, altLength)
  return (
    kind && {
      kind,
      refStart: d.refStart,
      refEnd: d.refEnd,
      altLength,
      delta: d.delta,
      ...(kind === 'duplication' ? { repeatedBp: d.repeatedBp } : {}),
      ...(kind === 'duplication' && d.copied ? { copied: d.copied } : {}),
    }
  )
}

/**
 * A deletion of 1 kb or more and an insertion or duplication within 10% of
 * its size, no further from it than twice that size, are R's length with the
 * sequence placed elsewhere: at AMY1, 191 haplotypes pair a 94,204 bp
 * deletion with a 94,142 bp insertion 18 kb away
 */
export function pairShifts(marks: Mark[]) {
  const used = new Set<Mark>()
  const shifts: Mark[] = []
  const deletions = marks
    .filter(m => m.kind === 'deletion' && -m.delta >= SHIFT_MIN_BP)
    .sort((a, b) => a.delta - b.delta)
  for (const d of deletions) {
    const size = -d.delta
    const gain = marks.find(
      m =>
        !used.has(m) &&
        (m.kind === 'insertion' || m.kind === 'duplication') &&
        Math.abs(m.delta - size) <= SHIFT_TOLERANCE * Math.max(m.delta, size) &&
        Math.max(m.refStart - d.refEnd, d.refStart - m.refEnd) <= 2 * size,
    )
    if (gain) {
      used.add(gain)
      used.add(d)
      const refStart = Math.min(d.refStart, gain.refStart)
      const refEnd = Math.max(d.refEnd, gain.refEnd)
      const delta = d.delta + gain.delta
      shifts.push({
        kind: 'shift',
        refStart,
        refEnd,
        altLength: refEnd - refStart + delta,
        delta,
        shiftedBp: size,
      })
    }
  }
  return [...marks.filter(m => !used.has(m)), ...shifts].sort(
    (a, b) => a.refStart - b.refStart,
  )
}

/**
 * A haplotype's marks of at least `floorBp`: its departures nearer than
 * `joinBp` on R joined, those touching the window kept, and shifts paired
 */
export function haplotypeMarks(
  departures: Departure[],
  {
    floorBp = STRUCTURAL_BP,
    joinBp = STRUCTURAL_BP,
    window,
  }: { floorBp?: number; joinBp?: number; window?: FormsOptions['window'] },
) {
  const marks: Mark[] = []
  for (const d of join(departures, joinBp)) {
    if (window && (d.refEnd < window.start || d.refStart > window.end)) {
      continue
    }
    const mark = markOf(d, floorBp)
    if (mark) {
      marks.push(mark)
    }
  }
  return pairShifts(marks)
}

const sumDelta = (marks: Mark[]) => marks.reduce((s, m) => s + m.delta, 0)

const covers = (h: HaplotypeDepartures, start: number, end: number) =>
  h.covered.some(([a, b]) => a <= start && b >= end)

const overlaps = (h: HaplotypeDepartures, start: number, end: number) =>
  h.covered.some(([a, b]) => a <= end && b >= start)

function siteSpans(marksBy: Mark[][]) {
  const spans = marksBy
    .flat()
    .map(m => [m.refStart, m.refEnd] as const)
    .sort((a, b) => a[0] - b[0])
  const sites: [number, number][] = []
  for (const [start, end] of spans) {
    const last = sites.at(-1)
    if (last && start <= last[1]) {
      last[1] = Math.max(last[1], end)
    } else {
      sites.push([start, end])
    }
  }
  return sites
}

function siteIndex(sites: [number, number][], m: Mark) {
  let lo = 0
  let hi = sites.length - 1
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1
    if (sites[mid]![0] <= m.refStart) {
      lo = mid
    } else {
      hi = mid - 1
    }
  }
  return lo
}

/**
 * Groups haplotypes by the structural changes their walks make over a
 * window. Each haplotype's own marks make its state at a site, so a change
 * nested in another counts once. A state with fewer than minCarriers folds
 * to the common state nearest in net change, and so does a form.
 */
export function structuralForms(
  input: HaplotypeDepartures[],
  opts: FormsOptions = {},
): StructuralForms {
  const exclude = opts.exclude ?? []
  const haplotypes = input.filter(
    h => !exclude.some(prefix => panSNMatchesPrefix(h.haplotype, prefix)),
  )
  const minCarriers = opts.minCarriers ?? minCarriersFor(haplotypes.length)
  const window = opts.window ?? extent(haplotypes)
  const floorBp =
    opts.floorBp === undefined || opts.floorBp === 'auto'
      ? autoFloor(window.end - window.start)
      : Math.max(STRUCTURAL_BP, opts.floorBp)
  const joinBp = opts.joinBp ?? STRUCTURAL_BP
  const marksBy = haplotypes.map(h =>
    haplotypeMarks(h.departures, { floorBp, joinBp, window: opts.window }),
  )
  const spans = siteSpans(marksBy)
  const marksAt = marksBy.map(marks => {
    const at = new Map<number, Mark[]>()
    for (const m of marks) {
      const i = siteIndex(spans, m)
      at.set(i, [...(at.get(i) ?? []), m])
    }
    return new Map(
      [...at].map(([site, ms]) => {
        const state = ms.map(stateKey).join(';')
        return [site, ms.map((m): SiteMark => ({ ...m, site, state }))]
      }),
    )
  })
  const raw = spans.map(([start, end], i) =>
    haplotypes.map((h, j) =>
      covers(h, Math.max(start, window.start), Math.min(end, window.end))
        ? (marksAt[j]!.get(i)?.[0]!.state ?? REFERENCE_STATE)
        : NO_CALL_STATE,
    ),
  )
  const inWindow = haplotypes.map(h => overlaps(h, window.start, window.end))
  const windowCount = inWindow.filter(Boolean).length
  const edge = raw.map(
    states =>
      2 * states.filter((s, j) => inWindow[j] && s !== NO_CALL_STATE).length <
      windowCount,
  )
  const isPlaced = haplotypes.map(
    (_, j) =>
      inWindow[j]! &&
      (edge.every(Boolean) ||
        raw.some((states, i) => !edge[i] && states[j] !== NO_CALL_STATE)),
  )
  const placed = haplotypes.flatMap((_, j) => (isPlaced[j] ? [j] : []))
  const rare = new Set<string>()
  const sites: Site[] = []
  const informative: number[] = []
  const folded: string[][] = []
  // the marks of every placed haplotype's state at a site that tells none
  // apart
  const shared: SiteMark[] = []
  spans.forEach(([refStart, refEnd], i) => {
    const states = raw[i]!
    const carriers = new Map<string, number[]>()
    for (const j of placed) {
      const list = carriers.get(states[j]!)
      if (list) {
        list.push(j)
      } else {
        carriers.set(states[j]!, [j])
      }
    }
    const ranked = [...carriers].sort(
      (a, b) => b[1].length - a[1].length || (a[0] < b[0] ? -1 : 1),
    )
    const siteStates = ranked.flatMap(([key, js]) => {
      if (key === NO_CALL_STATE) {
        return []
      }
      const byDelta = js
        .map(j => ({ j, marks: marksAt[j]!.get(i) ?? [] }))
        .map(c => ({ ...c, delta: sumDelta(c.marks) }))
        .sort((a, b) => a.delta - b.delta)
      const median = byDelta[Math.floor(byDelta.length / 2)]!
      return [
        {
          key,
          delta: median.delta,
          marks: median.marks,
          carriers: js.map(j => haplotypes[j]!.haplotype),
        },
      ]
    })
    const tells =
      !edge[i] && ranked.length > 1 && ranked[1]![1].length >= minCarriers
    sites.push({
      refStart,
      refEnd,
      edge: edge[i]!,
      shared: !edge[i] && !tells,
      covered: placed.length - (carriers.get(NO_CALL_STATE)?.length ?? 0),
      states: siteStates,
    })
    if (edge[i]) {
      return
    }
    const majority = ranked[0]?.[0] ?? REFERENCE_STATE
    const common = siteStates.filter(s => s.carriers.length >= minCarriers)
    const fold = states.map((state, j) => {
      const n = carriers.get(state)?.length ?? 0
      if (n >= minCarriers || !isPlaced[j]) {
        return state
      }
      if (state !== majority) {
        rare.add(haplotypes[j]!.haplotype)
      }
      if (state !== NO_CALL_STATE && majority === NO_CALL_STATE) {
        return RARER_STATE
      }
      if (state !== NO_CALL_STATE && common.length > 0) {
        const delta = sumDelta(marksAt[j]!.get(i) ?? [])
        return nearest(common, s => s.delta, delta).key
      }
      return majority
    })
    if (tells) {
      informative.push(i)
      folded.push(fold)
    } else {
      for (const m of siteStates.find(st => st.key === majority)?.marks ?? []) {
        shared.push(m)
      }
    }
  })
  const stateOf = (i: number, key: string) =>
    sites[informative[i]!]!.states.find(s => s.key === key)
  const byKey = new Map<string, { states: string[]; members: string[] }>()
  for (const j of placed) {
    const states = folded.map(fold => fold[j]!)
    const key = states.join('|')
    const form = byKey.get(key)
    if (form) {
      form.members.push(haplotypes[j]!.haplotype)
    } else {
      byKey.set(key, { states, members: [haplotypes[j]!.haplotype] })
    }
  }
  const forms = [...byKey].map(([key, { states, members }]): Form => {
    const marks = states
      .flatMap((state, i) => stateOf(i, state)?.marks ?? [])
      .concat(shared)
      .sort((a, b) => a.refStart - b.refStart)
    return { members, key, states, marks, net: sumDelta(marks) }
  })
  const common =
    opts.foldRareForms === false
      ? forms
      : forms.filter(f => f.members.length >= minCarriers)
  const kept = common.length > 0 ? common : forms
  if (common.length > 0 && common !== forms) {
    for (const form of sortForms(forms)) {
      if (form.members.length < minCarriers) {
        const into = nearest(common, f => f.net, form.net)
        for (const m of form.members) {
          into.members.push(m)
          rare.add(m)
        }
      }
    }
  }
  return {
    sites,
    informative,
    forms: sortForms(kept),
    placed: placed.length,
    unplaced: haplotypes.flatMap((h, j) => (isPlaced[j] ? [] : [h.haplotype])),
    rareCarriers: haplotypes
      .filter(h => rare.has(h.haplotype))
      .map(h => h.haplotype),
    abridged: haplotypes
      .filter((_, j) => marksBy[j]!.some(m => m.kind === 'net'))
      .map(h => h.haplotype),
    minCarriers,
    floorBp,
  }
}

function extent(haplotypes: HaplotypeDepartures[]) {
  let start = Infinity
  let end = -Infinity
  for (const h of haplotypes) {
    for (const [a, b] of h.covered) {
      start = Math.min(start, a)
      end = Math.max(end, b)
    }
  }
  return start <= end ? { start, end } : { start: 0, end: 0 }
}

// The first of `items` whose value is nearest `target`
function nearest<T>(items: T[], value: (item: T) => number, target: number) {
  let best: T = items[0]!
  for (const item of items) {
    if (Math.abs(value(item) - target) < Math.abs(value(best) - target)) {
      best = item
    }
  }
  return best
}

const sortForms = (forms: Form[]) =>
  forms.toSorted(
    (a, b) =>
      b.members.length - a.members.length ||
      Math.abs(a.net) - Math.abs(b.net) ||
      (a.key < b.key ? -1 : 1),
  )

export const bp = (n: number) =>
  n >= 1_000_000
    ? `${Number((n / 1_000_000).toPrecision(2))} Mb`
    : n >= 1000
      ? `${Number((n / 1000).toPrecision(2))} kb`
      : `${n} bp`

const bp3 = (n: number) =>
  n >= 1_000_000
    ? `${Number((n / 1_000_000).toPrecision(3))} Mb`
    : n >= 1000
      ? `${Number((n / 1000).toPrecision(3))} kb`
      : `${n} bp`

// A replacement far wider than its change says only how much longer or
// shorter the haplotype is there: AMY1's 176 kb site 2.8 kb longer
const DIFFUSE_MIN_SPAN = 10_000

const sizeOf = (m: Mark) => m.shiftedBp ?? Math.abs(m.delta)

/**
 * A mark in words: "18 kb deletion", "30.7 kb replaced by 109 kb", "49 kb
 * tandem duplication", "18 kb shifted, same length"
 */
export function changeName(m: Mark) {
  const span = m.refEnd - m.refStart
  const net = `${bp(Math.abs(m.delta))} ${m.delta > 0 ? 'longer' : 'shorter'} in ${bp(span)}`
  switch (m.kind) {
    case 'inversion':
      return `${bp(span)} inversion`
    case 'shift':
      return `${bp(m.shiftedBp!)} shifted, same length`
    case 'duplication':
      return `${bp(m.delta)} tandem duplication`
    case 'net':
      return net
    case 'replacement':
      return span > Math.max(2 * Math.abs(m.delta), DIFFUSE_MIN_SPAN)
        ? net
        : `${bp3(span)} replaced by ${bp3(m.altLength)}`
    default:
      return `${bp(Math.abs(m.delta))} ${m.kind}`
  }
}

const LISTED_CHANGES = 2

const sites = (n: number) => `${n} site${n === 1 ? '' : 's'}`

/**
 * How a form differs from R. Up to two size changes are listed, largest
 * first; more are counted with the largest named, since a list of sizes says
 * nothing.
 */
export function describeForm(form: { marks: Mark[]; states: string[] }) {
  const sized = form.marks
    .filter(m => m.kind !== 'inversion')
    .sort((a, b) => sizeOf(b) - sizeOf(a))
  const inversions = form.marks.length - sized.length
  const rarer = form.states.filter(s => s === RARER_STATE).length
  const uncalled = form.states.filter(s => s === NO_CALL_STATE).length
  const changes = [
    ...(sized.length > LISTED_CHANGES
      ? [`${sized.length} size changes, largest a ${changeName(sized[0]!)}`]
      : sized.map(changeName)),
    ...(inversions === 1
      ? [changeName(form.marks.find(m => m.kind === 'inversion')!)]
      : inversions > 1
        ? [`${inversions} inversions`]
        : []),
    ...(rarer > 0
      ? [
          rarer === 1
            ? 'a rarer change at 1 site'
            : `rarer changes at ${rarer} sites`,
        ]
      : []),
    ...(uncalled > 0 ? [`not aligned at ${sites(uncalled)}`] : []),
  ]
  return changes.length > 0 ? changes.join(', ') : 'as the reference'
}

/**
 * The state key `form` has at `site`, an index into `result.sites`: its key
 * there at an informative site, the state every form carries at a shared one,
 * and undefined at an edge site, which no form key holds
 */
export function formStateAt(
  result: Pick<StructuralForms, 'sites' | 'informative'>,
  form: Pick<Form, 'states' | 'marks'>,
  site: number,
) {
  const k = result.informative.indexOf(site)
  if (k >= 0) {
    return form.states[k]
  }
  return result.sites[site]?.shared
    ? (form.marks.find(m => m.site === site)?.state ?? REFERENCE_STATE)
    : undefined
}

const signed = (n: number) => `${n < 0 ? '\u2212' : '+'}${bp(Math.abs(n))}`

/**
 * A mark in a few characters: "+78 kb", "18 kb del", "49 kb dup", "2 kb
 * shift", "3 kb inv". An abridged walk gives its net change, "−94 kb", and so
 * does a replacement, unless one spanning 10 kb or less keeps most of its
 * length: "6 kb repl".
 */
export function markTag(m: Mark) {
  const span = m.refEnd - m.refStart
  switch (m.kind) {
    case 'insertion':
    case 'net':
      return signed(m.delta)
    case 'deletion':
      return `${bp(-m.delta)} del`
    case 'duplication':
      return `${bp(m.delta)} dup`
    case 'inversion':
      return `${bp(span)} inv`
    case 'shift':
      return `${bp(m.shiftedBp!)} shift`
    case 'replacement':
      return 2 * Math.abs(m.delta) < span && span <= DIFFUSE_MIN_SPAN
        ? `${bp(span)} repl`
        : signed(m.delta)
  }
}

const tagSize = (m: Mark) =>
  m.kind === 'inversion' ? m.refEnd - m.refStart : sizeOf(m)

/**
 * A form for a row label of about 20 characters: its largest change as
 * markTag gives it, then "+N" counting the rest describeForm names, as in
 * "18 kb del" or "+78 kb +1"; "reference" where it makes no change. Given
 * the form's result, the tag leaves out the marks every form shares, so a
 * form with no other change is "shared change": naming the change there
 * would read as if the other rows lacked it.
 */
export function formTag(
  form: { marks: (Mark | SiteMark)[]; states: string[] },
  result?: Pick<StructuralForms, 'sites'>,
) {
  const isShared = (m: Mark | SiteMark) =>
    'site' in m && result?.sites[m.site]?.shared === true
  const own = form.marks.filter(m => !isShared(m))
  const [largest] = own.toSorted((a, b) => tagSize(b) - tagSize(a))
  const shared = form.marks.length - own.length
  const rarer = form.states.filter(s => s === RARER_STATE).length
  const uncalled = form.states.filter(s => s === NO_CALL_STATE).length
  const head = largest
    ? markTag(largest)
    : rarer > 0
      ? 'rarer change'
      : uncalled > 0
        ? 'not aligned'
        : shared > 1
          ? 'shared changes'
          : shared === 1
            ? 'shared change'
            : 'reference'
  const rest = own.length + rarer + uncalled - 1
  return rest > 0 ? `${head} +${rest}` : head
}
