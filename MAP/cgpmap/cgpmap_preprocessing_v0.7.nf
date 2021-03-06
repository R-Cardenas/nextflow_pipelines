/*
 * create a channel for fastq pairs
 */

// Input Reads
//params.read1 = "/gpfs/afm/cg_pipelines/Datasets/Cholesteatoma batch 2/release_noclean/RawData/**/*{1,2}.fq.gz"
params.read1 = "$baseDir/input/*{1,2}.fq.gz"

read1_ch = Channel .fromFilePairs( params.read1 )
read1_ch.into { read2_ch; read3_ch }

params.csv = "$baseDir/bin/williams_batch2_info.csv"
csv_ch = Channel .fromPath( params.csv )


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
	file "${reads[0].simpleName}.fq.gz" into (read5_ch, read7_ch)
	file "${reads[1].simpleName}.fq.gz" into (read10_ch, read12_ch)
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

process fqtools{
	storeDir "$baseDir/output/cgpMAP/${read1.simpleName}/fqtools-config"
	input:
	file read1 from read7_ch
	file read2 from read12_ch
	output:
	file "*config.yaml" into yaml_ch
	script:
	"""
	fqtools -d header ${read1} | grep ":[C,A,T,G]*[+][C,A,T,G]" | head -1 > ${read1.simpleName}.txt
	fqtools -d header ${read2} | grep ":[C,A,T,G]*[+][C,A,T,G]" | head -1 > ${read2.simpleName}.txt

	python $baseDir/fastq2config_cgpmap.py --fq1 ${read1.simpleName}.txt --fq2 ${read1.simpleName}.txt \
	--n1 ${read1} --n2 ${read2}
	"""
}

process cgpMAP {
  storeDir "$baseDir/output/cgpMAP/${read1.simpleName}"
  input:
	val read1 from read5_ch
	val read2 from read10_ch
	val yaml from yaml_ch
  output:
  file "*.bam" into cgp_ch
  script:
  """
	echo 'fq1: ${read1} fq2: ${read2} bam_name: ${read1.simpleName}' >> cgpmap_samples.log
	name=\$(sed -n -e 's/^.*SM: //p' ${yaml})

  ds-cgpmap.pl  \
  -outdir $baseDir/output/cgpMAP/${read1.simpleName} \
  -r $cgpmap_genome \
  -i $cgpmap_index \
  -s \$name \
  -t 5 \
	-g $yaml \
  ${read1} ${read2}

	mv \$name.bam ${read1.simpleName}.bam 
  """
}


process sam_sort {
  storeDir "$baseDir/output/trim/aligned_sorted"
  input:
  file bam from cgp_ch
  output:
  file "${bam.simpleName}.sorted.bam" into dup_ch
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
  file bam from dup_ch
  output:
  file "${bam.simpleName}.rmd.bam" into (index1_ch, index_2ch, hs_ch, bam10_ch, bam11_ch)
	file "${bam.simpleName}.log"
  script:
  """
	mkdir -p tmp
  picard MarkDuplicates I=${bam} O=${bam.simpleName}.rmd.bam M=${bam.simpleName}.log TMP_DIR=tmp
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

process merge_lanes{
	storeDir "$baseDir/output/trim/merge_lanes"
	input:
	file bam from bam11_ch.collect()
	output:
	file "*_merge.bam" into rename2_ch
	script:
	"""
	module add python/anaconda/4.2/3.5
	$baseDir/bin/python merge_bam_lanes.py
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
	 merge bams - completed
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
