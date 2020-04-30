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
// link for adding fastp if the pipe works https://github.com/nextflow-io/nextflow/issues/682
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
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from cgp_ch
  output:
  file "${bam.simpleName}.sorted.bam" into sort_ch
  script:
  """
  picard SortSam I=${bam} O=${bam.simpleName}.sorted.bam SORT_ORDER=coordinate
  """
}

process picard_pcr_removal {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from sort_ch
  output:
  file "${bam.simpleName}.rmd.bam" into rename_ch
  script:
  """
  picard MarkDuplicates I=${bam} O=${bam.simpleName}.rmd.bam M=${bam.simpleName}.log
  """
}


process rename_bam {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from rename_ch
  output:
  file "${bam.simpleName}.rename.bam" into (index1_ch, index_2ch, hs_ch, bam10_ch)
  script:
  """
  picard AddOrReplaceReadGroups I=${bam} O=${bam.simpleName}.rename.bam \
  RGID=rename RGLB=${bam.simpleName} RGPL=illumina RGPU=unit1 RGSM=${bam.simpleName}
  """
}

process bam_index {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from index1_ch
  output:
  file "${bam}.bai"

  script:
  """
  samtools index ${bam}

  """
}

process collect_insert_size {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from index_2ch
  output:
  file "${bam.simpleName}_insert_size.txt"
  script:
  """
  picard CollectInsertSizeMetrics I=${bam} H=${bam.simpleName}_histogram.pdf \
  O=${bam.simpleName}_insert_size.txt M=0.5
  """
}

process hybrid_stats {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from hs_ch
  output:
  file "${bam.simpleName}_hs_metrics.txt"
  script:
  """
  picard CollectHsMetrics I=${bam} O=${bam.simpleName}_hs_metrics.txt \
  R=$genome_fasta \
  BAIT_INTERVALS=$bait_interval \
  TARGET_INTERVALS=$target_interval
  """
}

process alignment_stats{
  storeDir "$baseDir/output/trim/aligned_sorted"
	input:
	file bam from bam10_ch
	output:
	file "${bam.simpleName}_align_stats.txt"
	script:
	"""
	picard CollectAlignmentSummaryMetrics \
  R=$genome_fasta \
	I=${bam} \
	O=${bam.simpleName}_align_stats.txt
	"""
}

workflow.onComplete {
	// create log files and record output
	process finish{
		script:
		""" echo ' Pipeline cgpmap v0.3 completed
		Project: $projectname
		Time: ${nextflow.timestamp}

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
	 ' >> $baseDir/${projectname}_log.txt

	 mail -s "cgpMAP successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_log.txt
	 """
	}
}

workflow.onError {
	process finish_error{
		script:
		"""
		echo 'Pipeline cgpmap v0.3 FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}

		Error:
		${workflow.errorMessage}' >> $baseDir/${projectname}_error.txt

	  mail -s "cgpMAP successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_error.txt
	  """
	}
}
