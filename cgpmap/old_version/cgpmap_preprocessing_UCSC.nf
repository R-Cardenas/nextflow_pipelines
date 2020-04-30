/*
 * create a channel for fastq pairs
 */
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
         v0.1
         ===================================



         """
         .stripIndent()

process make_dir {
  input:
  tuple val(read2), file(reads) from read2_ch
  script:
  output:
  val "$baseDir/output/UCSC/cgpMAP/${read2}" into path_ch
  """
  mkdir -p $baseDir/output/UCSC/cgpMAP/${read2}
  """
}

process cgpMAP {
  storeDir "$baseDir/output/UCSC/cgpMAP/$read1"
  input:
  tuple val(read1), file(reads) from read3_ch
  val path from path_ch
  output:
  file "${read1}.bam" into cgp_ch

  script:
  """
  ds-cgpmap.pl  \
  -outdir ${path} \
  -r /var/spool/mail/hg38/UCSC/WholeGenome_hg38.tar.gz \
  -i /var/spool/mail/hg38/UCSC/BWAIndex_hg38.tar.gz \
  -s ${read1} \
  -t 5 \
  ${reads[0]} ${reads[1]}
  """
}

process sam_sort {
  storeDir "$baseDir/output/UCSC/aligned_sorted"
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
  storeDir "$baseDir/output/UCSC/aligned_sorted"
  input:
  file pair_read_3 from sort_ch
  output:
  file "${pair_read_3.simpleName}.rmd.bam" into (rename_ch, index66_ch)
  script:
  """
  picard MarkDuplicates I=${pair_read_3} O=${pair_read_3.simpleName}.rmd.bam M=${pair_read_3.simpleName}.log
  """
}



process bam_index {
  storeDir "$baseDir/output/UCSC/aligned_sorted"
  input:
  file pair_read_8 from index66_ch
  output:
  file "${pair_read_8.simpleName}.bai"

  script:
  """
  samtools index ${pair_read_8} ${pair_read_8.simpleName}.bai

  """
}

process rename_bam {
  storeDir "$baseDir/output/UCSC/aligned_sorted"
  input:
  file pair_read_7 from rename_ch
  output:
  file "${pair_read_7.simpleName}.rename.bam" into (index_2ch, hs_ch)
  script:
  """
  picard AddOrReplaceReadGroups I=${pair_read_7} O=${pair_read_7.simpleName}.rename.bam \
  RGID=rename RGLB=${pair_read_7.simpleName} RGPL=illumina RGPU=unit1 RGSM=${pair_read_7.simpleName}
  """
}

process collect_insert_size {
  storeDir "$baseDir/output/UCSC/aligned_sorted"
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
  storeDir "$baseDir/output/UCSC/aligned_sorted"
  input:
  file bam11 from hs_ch
  output:
  file "${bam11.simpleName}_hs_metrics.txt"
  script:
  """
  picard CollectHsMetrics I=${bam11} O=${bam11.simpleName}_hs_metrics.txt \
  R=/var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  BAIT_INTERVALS=/var/spool/mail/bait_capture_files/SeqCap_EZ_MedExome/MedExome_hg38_capture_targets_UCSC.interval \
  TARGET_INTERVALS=/var/spool/mail/bait_capture_files/SeqCap_EZ_MedExome/MedExome_hg38_empirical_targets_UCSC.interval
  """
}
