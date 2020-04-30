/*
 * create a channel for fastq pairs
 */

// Input Reads
params.read1 = "/gpfs/afm/cg_pipelines/Datasets/Cholesteatoma batch 2/release_noclean/RawData/**/*{1,2}.fq.gz"
read1_ch = Channel .fromFilePairs( params.read1 )


println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         Mapping and Bam Processing
         v0.3
         ===================================



         """
         .stripIndent()


// link for adding fastp if the pipe works https://github.com/nextflow-io/nextflow/issues/682
process trimmomatic{
  storeDir "$baseDir/output/cgpMAP/trimmomatic"
	input:
	tuple val(read2), file(reads) from read1_ch
	output:
	file "${reads[0].simpleName}_1.fq.gz" into read4_ch
	file "${reads[1].simpleName}_2.fq.gz" into (read5_ch, read6_ch, read7_ch)
	script:
	"""
	module add trimgalore/0.4.3
	trim_galore --paired ${reads[0]} ${reads[1]} \
	--fastqc -a $adapter1 \
	-a2 $adapter2

	mv ${reads[0].simpleName}_val_1.fq.gz ${reads[0].simpleName}_1.fq.gz
	mv ${reads[1].simpleName}_val_2.fq.gz ${reads[1].simpleName}_2.fq.gz
	"""
}

process make_dir {
  input:
  val read2 from read6_ch
  output:
  val "$baseDir/output/cgpMAP/${read2}" into path_ch
  script:
  """
  mkdir -p $baseDir/output/cgpMAP/${read2}
  """
}

process cgpMAP {
  storeDir "$baseDir/output/cgpMAP/$x"
  input:
	file read1 from read4_ch
	file read2 from read7_ch
	val x from read5_ch
	val path from path_ch
  output:
  file "${read1.simpleName}.bam" into cgp_ch
  script:
  """
  ds-cgpmap.pl  \
  -outdir ${path} \
  -r $cgpmap_genome \
  -i $cgpmap_index \
  -s ${read1.simpleName} \
  -t 5 \
  ${read1} ${read2}
  """
}

// merge bam here from different lanes

process sam_sort {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from cgp_ch
  output:
  file "${bam.simpleName}.sorted.bam" into sort_ch
  script:
  """
	mkdir $baseDir/tmp
  picard SortSam I=${bam} O=${bam.simpleName}.sorted.bam SORT_ORDER=coordinate TMP_DIR=$baseDir/tmp
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
  picard MarkDuplicates I=${bam} O=${bam.simpleName}.rmd.bam M=${bam.simpleName}.log TMP_DIR=$baseDir/tmp
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
  RGID=rename RGLB=${bam.simpleName} RGPL=illumina RGPU=unit1 RGSM=${bam.simpleName} TMP_DIR=$baseDir/tmp
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
  picard BuildBamIndex \
	I=${bam} \
	TMP_DIR=$baseDir/tmp

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
  O=${bam.simpleName}_insert_size.txt M=0.5 TMP_DIR=$baseDir/tmp
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
  TARGET_INTERVALS=$target_interval \
	TMP_DIR=$baseDir/tmp
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
	O=${bam.simpleName}_align_stats.txt \
	TMP_DIR=$baseDir/tmp
	"""
}

workflow.onComplete {
	// create log files and record output
	process finish{
		script:
		""" echo ' Pipeline cgpmap v0.3 completed
		Project: $projectname
		Time: ${nextflow.timestamp}

	 trimmomatic - completed
	 (adapter1: $adapter1
		adapter2: $adapter2)
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
