/*
 * create a channel for fastq pairs
 */

// Input Reads
params.read1 = "/gpfs/afm/cg_pipelines/Datasets/Cholesteatoma batch 2/release_noclean/RawData/**/*{1,2}.fq.gz"
read1_ch = Channel .fromFilePairs( params.read1 )
read1_ch.into { read2_ch; read3_ch }

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
process trim_galore{
  storeDir "$baseDir/output/cgpMAP/trimmomatic"
	input:
	tuple val(read2), file(reads) from read2_ch
	output:
	file "${reads[0].simpleName}.fq.gz" into read5_ch
	file "${reads[1].simpleName}.fq.gz" into read10_ch
	script:
	"""
	module add trimgalore/0.4.3
	trim_galore --paired ${reads[0]} ${reads[1]} \
	--fastqc -a $adapter1 \
	-a2 $adapter2

	mv ${reads[0].simpleName}_val_1.fq.gz ${reads[0].simpleName}.fq.gz
	mv ${reads[1].simpleName}_val_2.fq.gz ${reads[1].simpleName}.fq.gz
	"""
}


process cgpMAP {
  storeDir "$baseDir/output/cgpMAP/${read1.simpleName}"
  input:
	val read1 from read5_ch
	val read2 from read10_ch
  output:
  file "${read1.simpleName}.bam" into cgp_ch
  script:
  """
	echo 'fq1: ${read1} fq2: ${read2} bam_name: ${read1.simpleName}'

  ds-cgpmap.pl  \
  -outdir $baseDir/output/cgpMAP/${read1.simpleName} \
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
	mkdir -p tmp
  picard SortSam I=${bam} O=${bam.simpleName}.sorted.bam SORT_ORDER=coordinate TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process picard_pcr_removal {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from sort_ch
  output:
  file "${bam.simpleName}.rmd.bam" into rename_ch
	file "${bam.simpleName}.log"
  script:
  """
	mkdir -p tmp
  picard MarkDuplicates I=${bam} O=${bam.simpleName}.rmd.bam M=${bam.simpleName}.log TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
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
	mkdir tmp
  picard AddOrReplaceReadGroups I=${bam} O=${bam.simpleName}.rename.bam \
  RGID=rename RGLB=${bam.simpleName} RGPL=illumina RGPU=unit1 RGSM=${bam.simpleName} TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
	"""
}

process bam_index {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from index1_ch
  output:
  file "${bam.simpleName}.rename.bai"

  script:
  """
	mkdir -p tmp
  picard BuildBamIndex \
	I=${bam} \
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process collect_insert_size {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from index_2ch
  output:
  file "${bam.simpleName}_insert_size.txt"
	file "${bam.simpleName}_insert_size.txt"
  script:
  """
	mkdir -p tmp
  picard CollectInsertSizeMetrics I=${bam} H=${bam.simpleName}_histogram.pdf \
  O=${bam.simpleName}_insert_size.txt M=0.5 TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
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
	mkdir -p tmp
  picard CollectHsMetrics I=${bam} O=${bam.simpleName}_hs_metrics.txt \
  R=$genome_fasta \
  BAIT_INTERVALS=$bait_interval \
  TARGET_INTERVALS=$target_interval \
	TMP_DIR=tmp
	rm -fr tmp
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
	mkdir -p tmp
	picard CollectAlignmentSummaryMetrics \
  R=$genome_fasta \
	I=${bam} \
	O=${bam.simpleName}_align_stats.txt \
	TMP_DIR=tmp
	rm -fr tmp
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
