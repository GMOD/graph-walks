import type { GraphTables } from './util.ts'

// The GFA a cut's tables stand for, as @jbrowse/bandage-core writes it
export function graphTablesGFA({ nodes, links, walks }: GraphTables) {
  const { names } = nodes
  const lines = ['H\tVN:Z:1.1']
  for (let i = 0; i < nodes.lengths.length; i++) {
    const ref = nodes.refs[i]!
    const stable =
      ref < 0
        ? ''
        : `\tSN:Z:${nodes.refNames[ref]}\tSO:i:${nodes.starts[i]}\tSR:i:${nodes.ranks[i]}`
    lines.push(`S\t${names[i]}\t*\tLN:i:${nodes.lengths[i]}${stable}`)
  }
  for (let i = 0; i < links.from.length; i++) {
    const strands = links.strands[i]!
    lines.push(
      `L\t${names[links.from[i]!]}\t${strands & 1 ? '-' : '+'}\t${names[links.to[i]!]}\t${strands & 2 ? '-' : '+'}\t0M`,
    )
  }
  for (let p = 0; p < walks.names.length; p++) {
    const [sample, haplotype, ...contig] = walks.names[p]!.split('#')
    let body = ''
    for (let s = walks.offsets[p]!; s < walks.offsets[p + 1]!; s++) {
      body += (walks.reversed[s] ? '<' : '>') + names[walks.steps[s]!]
    }
    lines.push(
      `W\t${sample}\t${haplotype}\t${contig.join('#')}\t${walks.starts[p]}\t${walks.ends[p]}\t${body}`,
    )
  }
  return lines.join('\n')
}
