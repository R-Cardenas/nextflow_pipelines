/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/input/*.bam"
bam_ch = Channel .fromPath( params.bam )

println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         Picard SamToFastQ
				 v0.1
         ===================================



         """
         .stripIndent()

process SamToFastQ{
  storeDir "$baseDir/output"
  input:
  file bam from bam_ch
  output:
  file "${bam.simpleName}_1.fq.gz"
  file "${bam.simpleName}_2.fq.gz"
  script:
  """
	samtools bam2fq ${bam} > ${bam.simpleName}.fastq

	cat ${bam.simpleName}.fastq  | grep '^@.*/1\$' -A 3 --no-group-separator > ${bam.simpleName}_1.fq
	cat ${bam.simpleName}.fastq  | grep '^@.*/2\$' -A 3 --no-group-separator > ${bam.simpleName}_2.fq

	gzip ${bam.simpleName}_1.fq
	gzip ${bam.simpleName}_2.fq

  """
}
