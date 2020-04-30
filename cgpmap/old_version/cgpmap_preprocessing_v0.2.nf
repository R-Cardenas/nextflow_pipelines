/*
 * create a channel for fastq pairs
 */

// Input Reads
params.read1 = "$baseDir/input/*{1,2}.fq.gz"
read1_ch = Channel .fromFilePairs( params.read1 )
read1_ch.into { read2_ch; read3_ch }


println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         Mapping and Bam Processing
         v0.2
         ===================================



         """
         .stripIndent()

process make_dir {
  input:
  tuple val(read2), file(reads) from read2_ch
  output:
  val "$baseDir/output/cgpMAP/${read2}" into path_ch
  script:
  """
  mkdir -p $baseDir/output/cgpMAP/${read2}
  """
}

process cgpMAP {
  storeDir "$baseDir/output/cgpMAP/$read2"
  input:
  tuple val(read2), file(reads) from read3_ch
  val path from path_ch

  output:
  file "${read2}.bam" into cgp_ch

  script:
  """
  ds-cgpmap.pl  \
  -outdir ${path} \
  -r $cgpmap_genome \
  -i $cgpmap_index \
  -s ${read2} \
  -t 5 \
  ${reads[0]} ${reads[1]}
  """
}

process sam_sort {
  storeDir "$baseDir/output/hg38_decoy/aligned_sorted"
  input:
  file pair_read_2 from cgp_ch
  output:
  file "${pair_read_2.simpleName}.sorted.bam" into sort_ch
  script:
  """
  picard SortSam I=${pair_read_2} O=${pair_read_2.simpleName}.sorted.bam SORT_ORDER=coordinate
  """
}

process picard_pcr_removal {
  storeDir "$baseDir/output/hg38_decoy/aligned_sorted"
  input:
  file pair_read_3 from sort_ch
  output:
  file "${pair_read_3.simpleName}.rmd.bam" into rename_ch
  script:
  """
  picard MarkDuplicates I=${pair_read_3} O=${pair_read_3.simpleName}.rmd.bam M=${pair_read_3.simpleName}.log
  """
}


process rename_bam {
  storeDir "$baseDir/output/hg38_decoy/aligned_sorted"
  input:
  file pair_read_7 from rename_ch
  output:
  file "${pair_read_7.simpleName}.rename.bam" into (index66_ch, index_2ch, hs_ch)
  script:
  """
  picard AddOrReplaceReadGroups I=${pair_read_7} O=${pair_read_7.simpleName}.rename.bam \
  RGID=rename RGLB=${pair_read_7.simpleName} RGPL=illumina RGPU=unit1 RGSM=${pair_read_7.simpleName}
  """
}

process bam_index {
  storeDir "$baseDir/output/hg38_decoy/aligned_sorted"
  input:
  file pair_read_8 from index66_ch
  output:
  file "${pair_read_8}.bai"

  script:
  """
  samtools index ${pair_read_8}

  """
}

process collect_insert_size {
  storeDir "$baseDir/output/hg38_decoy/aligned_sorted"
  input:
  file bam10 from index_2ch
  output:
  file "${bam10.simpleName}_insert_size.txt"
  script:
  """
  picard CollectInsertSizeMetrics I=${bam10} H=${bam10.simpleName}_histogram.pdf \
  O=${bam10.simpleName}_insert_size.txt M=0.5
  """
}

process hybrid_stats {
  storeDir "$baseDir/output/hg38_decoy/aligned_sorted"
  input:
  file bam11 from hs_ch
  output:
  file "${bam11.simpleName}_hs_metrics.txt"
  script:
  """
  picard CollectHsMetrics I=${bam11} O=${bam11.simpleName}_hs_metrics.txt \
  R=$genome_fasta \
  BAIT_INTERVALS=$bait_interval \
  TARGET_INTERVALS=$target_interval
  """
}

workflow.onComplete {
	// create log files and record output
	log = """ echo ' Pipeline cgpmap v0.2 completed

	cgpmap - completed
	(reference = $cgpmap_genome)
	(index = $cgpmap_index)
	samsort - completed
	picard pcr removal - completed
	picard rename bam - completed
	samtools bam index - completed
	picard collect insert size - completed
	picard collect HS stats - completed
	(target_intervals = $target_interval)
	(bait_intervals = $target_interval)
	' >> log.txt
	"""
	println log.execute().text

	// Send email user
	sendMail( to: 'aft19qdu@uea.ac.uk',
          subject: 'nextflow cgpmap pipeline',
          body: 'Dear User, cgpmap has now completed successfully. Please find log file attached.',
          attach: 'log.txt' )

}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"

		sendMail( to: 'aft19qdu@uea.ac.uk',
	          subject: 'nextflow cgpmap pipeline - FAIL',
	          body: 'Dear User, cgpmap has failed. Below is the error message: \n ${workflow.errorMessage}',
	          attach: 'log.txt' )
}
