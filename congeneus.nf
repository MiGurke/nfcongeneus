#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow consensus creater
         ===================================
         Bam directory           : ${params.bamdir}
         Reference annotation    : ${params.anno}
         Reference file          : ${params.ref}
         Chromosome list         : ${params.chr}
         Window size             : ${params.win_size}
         Feature                 : ${params.feat}
         Min Depth               : ${params.mindepth}
         Output directory        : ${params.outdir}
         """
         .stripIndent()

 // holds the individual bam files
 Channel
   .fromPath("${params.bamdir}*/*.bam")
   .set{bam_ch}

if (params.feat != null) {
 process GetFeat {
   output:
   stdout into chunk_list
   stdout into chunk_list2

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
    stdout into chunk_list2

    script:
    """
  #!/usr/bin/env python
  from Bio import SeqIO

  chr = "$chr"
  ref = SeqIO.parse("$params.ref","fasta")
  piece = $params.win_size

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
                      start = piece * i + 1
                      end = len(rec.seq)
                  else:
                      start = piece * (i-1) + 1
                      end = piece * i
                  line = "$chr"+":"+str(start)+"-"+str(end)
                  print(line)
    """
  }
}

chunk_ch = chunk_list.splitText()
chunk_ch2 = chunk_list2.splitText()

process IndexBAM {
  input:
  file(bam) from bam_ch

  output:
  tuple file(bam), file('*.bai') into bambai_ch

  script:
  """
  samtools index $bam
  """
}

process GetBams {

  input:
  tuple file(bam), file(bai) from bambai_ch
  each chunk from chunk_ch

  output:
  tuple file('*.bam'), val(chunkkkk) into featBam_ch

  script:
  chunkk = chunk.trim()
  chunkkk = chunkk.replaceAll(':','_')
  chunkkkk = chunkkk.replaceAll('\\.','_')
  """
  samtools view -h -X $bam $bai ${chunkk} > ${bam.baseName}_${chunkkkk}.bam
  """
}

process CreateConsensus {

  input:
  tuple file(bam), val(chunk) from featBam_ch

  output:
  tuple val(chunk), stdout into fasta_ch

  script:
  """
  angsd -doFasta 1 -doCounts 1 -setMinDepth ${params.mindepth} -i $bam
  zcat angsdput.fa.gz | sed 's/>.*/>${bam.baseName}/'
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
  echo "${fasta}" | sed -e 's/,//g' -e 's/, *//g' -e 's/ //g' -e 's/\\[//g' -e 's/\\]//g' > ${chunk}.fasta
  """

}
