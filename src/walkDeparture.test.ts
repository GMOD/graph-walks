import {
  adjacentSteps,
  departureKind,
  departureRowKey,
  departureStretch,
  parseDepartureRow,
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
    'GRCh38#0#chr1~departure\t19\t81\tdeparture\tHG002#1#ctg\t20\t30\t12\t20\t+\t80\t-'.split(
      '\t',
    )
  expect(parseDepartureRow(cols)).toEqual({
    name: 'HG002#1#ctg',
    hapStart: 20,
    ownBp: 30,
    missingBp: 12,
    leavePos: 20,
    leaveReversed: false,
    rejoinPos: 80,
    rejoinReversed: true,
    fromSteps: false,
  })
  expect(departureRowKey(cols)).toBe(
    'departure\tHG002#1#ctg\t20\t30\t12\t20\t+\t80\t-',
  )
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
