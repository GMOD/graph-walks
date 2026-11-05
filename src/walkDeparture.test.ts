import {
  adjacentSteps,
  departureKind,
  departureRowKey,
  departureStretch,
  parseDepartureRow,
  placementChoice,
  squeezer,
} from './walkDeparture.ts'

const ends = (
  leavePos: number,
  leave: '+' | '-',
  rejoinPos: number,
  rejoin: '+' | '-',
  ownBp = 0,
) => ({
  ownBp,
  leavePos,
  leaveReversed: leave === '-',
  rejoinPos,
  rejoinReversed: rejoin === '-',
})

test('a record row parses to its ends, own sequence and contig position', () => {
  const cols =
    'GRCh38#0#chr1~departure\t19\t81\tdeparture\tHG002#1#ctg\t20\t30\t12\t20\t+\t1\t80\t-\t3'.split(
      '\t',
    )
  expect(parseDepartureRow(cols)).toEqual({
    name: 'HG002#1#ctg',
    hapStart: 20,
    ownBp: 30,
    missingBp: 12,
    leavePos: 20,
    leaveReversed: false,
    leavePasses: 1,
    rejoinPos: 80,
    rejoinReversed: true,
    rejoinPasses: 3,
    fromSteps: false,
  })
  expect(placementChoice(parseDepartureRow(cols))).toBe(true)
  expect(departureRowKey(cols)).toBe(
    'departure\tHG002#1#ctg\t20\t30\t12\t20\t+\t1\t80\t-\t3',
  )
})

test('a record row without pass counts is refused, not misread', () => {
  const cols =
    'GRCh38#0#chr1~departure\t19\t81\tdeparture\tHG002#1#ctg\t20\t30\t12\t20\t+\t80\t-'.split(
      '\t',
    )
  expect(() => parseDepartureRow(cols)).toThrow(/not a walk format 6 row/)
})

test('a departure is a placement choice where either end lies on a node the reference passes more than once', () => {
  expect(placementChoice({ leavePasses: 1, rejoinPasses: 1 })).toBe(false)
  expect(placementChoice({ leavePasses: 2, rejoinPasses: 1 })).toBe(true)
  expect(placementChoice({ leavePasses: 1, rejoinPasses: 2 })).toBe(true)
})

test('a departure is a strand switch, a jump back, own sequence between adjacent steps or a jump forward', () => {
  expect(departureKind(ends(20, '+', 80, '-'))).toBe('switch')
  expect(departureKind(ends(80, '+', 20, '+'))).toBe('back')
  expect(departureKind(ends(20, '-', 80, '-'))).toBe('back')
  expect(departureKind(ends(20, '+', 20, '+', 30))).toBe('adjacent')
  expect(departureKind(ends(20, '+', 80, '+'))).toBe('forward')
  expect(departureKind(ends(80, '-', 20, '-'))).toBe('forward')
})

test('steps either side of a hole in the reference are adjacent', () => {
  const squeeze = squeezer(
    [
      { refName: 'R#0#chr', start: 12, end: 20 },
      { refName: 'R#0#chr2', start: 0, end: 50 },
    ],
    'R#0#chr',
  )
  expect([squeeze(5), squeeze(12), squeeze(20), squeeze(32)]).toEqual([
    5, 12, 12, 24,
  ])
  expect(adjacentSteps(ends(12, '+', 20, '+'), squeeze)).toBe(true)
  expect(adjacentSteps(ends(12, '+', 20, '+'), bp => bp)).toBe(false)
  expect(adjacentSteps(ends(12, '+', 20, '-'), squeeze)).toBe(false)
  expect(adjacentSteps(ends(12, '+', 20, '+', 8), squeeze)).toBe(false)
  expect(departureKind(ends(12, '+', 20, '+', 8), squeeze)).toBe('adjacent')
  expect(squeezer([], 'R#0#chr')(100)).toBe(100)
})

test('a departure covers the bases beside its breakpoints, as its record is filed', () => {
  expect(departureStretch(ends(20, '+', 80, '+'))).toEqual([19, 81])
  expect(departureStretch(ends(70, '-', 40, '-'))).toEqual([39, 71])
  expect(departureStretch(ends(70, '-', 50, '-'))).toEqual([49, 71])
  expect(departureStretch(ends(20, '+', 20, '+', 6))).toEqual([19, 21])
  expect(departureStretch(ends(0, '-', 0, '+'))).toEqual([0, 1])
})
