
/*
 * create a channel for fastq pairs
 */
params.bam = "$baseDir/output/aligned_sorted/*.rename.bam"
params.csv = "$baseDir/bin/williams_batch2_info.csv"
bam_ch = Channel .fromPath( params.bam )
csv_ch = Channel .fromPath( params.csv )

process merge_lanes{
	storeDir "$baseDir/output/trim/merge_lanes"
	input:
	file bam from bam_ch.collect()
  file csv from csv_ch
	output:
	file "*_merge.bam" into rename2_ch
	script:
	"""
	module add python/anaconda/4.2/3.5
	python $baseDir/bin/merge_bam_lanes.py
	"""
}

process rename_bam2 {
  storeDir "$baseDir/output/trim/merge_lanes"
  input:
  file bam from rename2_ch.flatMap()
  output:
  file "${bam.simpleName}.rename.bam"
  script:
  """
	mkdir tmp
  picard AddOrReplaceReadGroups I=${bam} O=${bam.simpleName}.rename.bam \
  RGID=rename RGLB=${bam.simpleName} RGPL=illumina RGPU=unit1 RGSM=${bam.simpleName} TMP_DIR=tmp
	rm -fr tmp
	"""
}
