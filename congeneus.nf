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
  skip = $params.skip

  for rec in ref:
      id = rec.id
      if id == chr:
          nparts = (len(rec.seq)/(piece+skip)) + 1
          if nparts < 1:
              line = "$chr"
              print(line)
          else:
              for i in range(0,int(nparts)):
                  if i == 0:
                      start = 0
                      end = piece
                  elif i == int(nparts) - 1:
                      start = (piece + skip) * i
                      if (start + piece) > len(rec.seq):
                        end = len(rec.seq)
                      else:
                        end = start + piece
                  else:
                      start = (piece + skip) * i
                      end = start + piece

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

process CreateConsensus {
  
  input:
  tuple file(bam), file(bai) from bambai_ch
  each chunk from chunk_ch

  output:
  tuple val(chunk), file('*.fasta') into fasta_ch

  script:
  chunkk = chunk.trim()
  chunkkk = chunkk.replaceAll(':','_')
  chunkkkk = chunkkk.replaceAll('\\.','_')
  """
  start=\$(echo $chunkk | sed 's/.*://' | sed 's/-.*//')
  end=\$(echo $chunkk | sed 's/.*-//')
  angsd -doFasta 1 -doCounts 1 -r ${chunkk} -setMinDepth ${params.mindepth} -i ${bam}
  python3 ${projectDir}/bin/red_fasta.py -f angsdput.fa.gz -s \$start -e \$end -b ${bam.baseName} -c ${chunkk} > ${bam.baseName}_${chunkkkk}.fasta
  """
}

mfasta_ch = fasta_ch.groupTuple()

process WriteFasta {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(chunk), file(fasta) from mfasta_ch

  output:
  file('*')

  script:
  """
  cat ${fasta} | sed -e 's/,//g' -e 's/, *//g' -e 's/ //g' -e 's/\\[//g' -e 's/\\]//g' > ${chunk.trim()}.fasta
  """

}
