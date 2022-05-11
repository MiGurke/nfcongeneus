#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow raxml phylogenies
         ===================================
         Fasta directory           : ${params.fastadir}
         Output directory          : ${params.outdir}
         ML searches               : ${params.MLsearches}
         """
         .stripIndent()

 // holds the individual fasta files
 Channel
   .fromPath("${params.fastadir}*.fa")
   .set{fasta_ch}

process raxml {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(fasta) from fasta_ch

  output:
  file('*')

  script:
  """
  raxmlHPC-PTHREADS -m GTRGAMMA -T ${task.cpus} -p 12345 -s $fasta -N ${params.MLsearches} -n  ${fasta.baseName}
  """
}
