#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow consensus creater
         ===================================
         Bam directory           : ${params.bamdir}
         Reference annotation    : ${params.anno}
         Reference file          : ${params.ref}
         Chromosome list         : ${params.chr}
         Window size             : ${params.win_size}
         Region skip             : ${params.skip}
         Feature                 : ${params.feat}
         Min Depth               : ${params.mindepth}
         Output directory        : ${params.outdir}
         """
         .stripIndent()

 // holds the individual bam files and their index files
 Channel
    .fromFilePairs("${params.bamdir}*/*{.bam,.bam.bai}")
    .set{bambai_ch}


if (params.feat != null) {
 process GetFeat {

   output:
   stdout into chunk_list

   script:
   """
   awk '\$3 == "${params.feat}" {print \$1":"\$4"-"\$5}' ${params.anno}
   """
 }
}

if (params.win_size != null) {

  chr_list = file(params.chr).readLines()

  process GetWin {

    input:
    val(chr) from chr_list

    output:
    stdout into chunk_list

    script:
    """
  #!/usr/bin/env python
  from Bio import SeqIO

  chr = "$chr"
  ref = SeqIO.parse("$params.ref","fasta")
  piece = $params.win_size
  skip = $params.skip

  for rec in ref:
      id = rec.id
      if id == chr:
          nparts = len(rec.seq)/piece
          if nparts < 1:
              line = "$chr"
              print(line)
          else:
              for i in range(1,int(nparts) + 1):
                  if i == 1:
                      start = 0
                      end = piece * i
                  elif i == int(nparts):
                      start = (piece + skip) * i + 1
                      end = len(rec.seq)
                  else:
                      start = (piece + skip) * (i-1) + 1
                      end = start + piece
                  line = "$chr"+":"+str(start)+"-"+str(end)
                  print(line)
    """
  }
}

chunk_ch = chunk_list.splitText()


process GetBams {

  input:
  set val(sample_id), file(bam), file(bai) from bambai_ch
  each chunk from chunk_ch

  output:
  tuple file('*.bam'), val(chunkkkk) into featBam_ch

  script:
  chunkk = chunk.trim()
  chunkkk = chunkk.replaceAll(':','_')
  chunkkkk = chunkkk.replaceAll('\\.','_')
  """
  samtools view -h -X $bam $bai ${chunkk} > ${sample_id}_${chunkkkk}.bam
  """
}

process CreateConsensus {

  input:
  tuple file(bam), val(chunk) from featBam_ch

  output:
  tuple val(chunk), file('consensus.fasta') into fasta_ch

  script:
  """
  angsd -doFasta 1 -doCounts 1 -setMinDepth ${params.mindepth} -i $bam
  zcat angsdput.fa.gz | sed 's/>.*/>${bam.baseName}/' > consensus.fasta
  """
}

mfasta_ch = fasta_ch.groupTuple()

process WriteFasta {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(chunk), val(fasta) from mfasta_ch

  output:
  file('*')

  script:
  """
  cat ${fasta} | sed -e 's/,//g' -e 's/, *//g' -e 's/ //g' -e 's/\\[//g' -e 's/\\]//g' > ${chunk}.fasta
  """

}
