/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "/gpfs/afm/cg_pipelines/Datasets/Cholesteatoma/Barbara_Jennings_UEA_BJ_ENQ-1485_A_02_ext_Analysis/bowtie2_sorted_bam_files/*.bam"
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
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output"
	stageInMode "copy"
  input:
  file bam from bam_ch
  output:
  file "${bam.simpleName}_1.fq.gz"
  file "${bam.simpleName}_2.fq.gz"
  script:
  """
	picard SamToFastq I=${bam} \
	F=${bam.simpleName}_1.fq F2=${bam.simpleName}_2.fq

	gzip *.fq
  """
}
