/*
 * create a channel for fastq pairs
 */
 params.read1 = "$baseDir/input/*1.fq.gz"
 params.read2 = "$baseDir/input/*2.fq.gz"

read1_ch = Channel .fromPath( params.read1 )
read2_ch = Channel .fromPath( params.read2 )


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

process FastP{
  storeDir "$baseDir/output/FastP"
  input:
  file read1 from read1_ch
  file read2 from read2_ch
  output:
  file "${read1.simpleName}_fp.fq.gz" into fastp1_ch
  file "${read2.simpleName}_fp.fq.gz" into fastp2_ch
  file "${read1.simpleName}_fastp.json"
  file "${read1.simpleName}_fastp.html"
  script:
  """
  fastp -i ${read1} -I ${read2} -o ${read1.simpleName}_fp.fq.gz -O ${read2.simpleName}_fp.fq.gz \
  --json ${read1.simpleName}_fastp.json \
  --html ${read1.simpleName}_fastp.html
  """
}

process cgpMAP {
  storeDir "$baseDir/output/cgpMAP"
  input:
  file read1 from fastp1_ch
  file read2 from fastp2_ch

  output:
  file "${read1.simpleName}.bam" into cgp_ch

  script:
  """
  ds-cgpmap.pl  \
  -r /var/spool/mail/hg38/UCSC/WholeGenome_hg38.tar.gz \
  -i /var/spool/mail/hg38/UCSC/BWAIndex_hg38.tar.gz \
  -s ${read1.simpleName}  \
  -t 15 \
  -outdir $baseDir/output/cgpMAP \
  ${read1} ${read2}
  """
}

process sam_sort {
  scratch true
  input:
  file pair_read_2 from cgp_ch
  output:
  file "${pair_read_2.simpleName}.sorted.bam" into sort_ch
  script:
  """
  samtools sort $pair_read_2 -O BAM -o ${pair_read_2.simpleName}.sorted.bam
  """
}


process picard_pcr_removal {
  scratch true
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
  storeDir "$baseDir/output/aligned_sorted"
  input:
  file pair_read_7 from rename_ch
  output:
  file "${pair_read_7.simpleName}_processed.bam" into bam_ch
  script:
  """
  picard AddOrReplaceReadGroups I=$pair_read_7 O=${pair_read_7.simpleName}_processed.bam \
  RGID=rename RGLB=${pair_read_7.simpleName} RGPL=illumina RGPU=unit1 RGSM=${pair_read_7.simpleName}
  """
}

process bam_index {
  storeDir "$baseDir/output/aligned_sorted"
  input:
  file pair_read_8 from bam_ch
  output:
  file "${pair_read_8.simpleName}.bai" into index_ch
  file "${pair_read_8.simpleName}_hs_metrics.txt"
  file "${pair_read_8.simpleName}_insert_size.txt"
  script:
  """
  samtools index ${pair_read_8} ${pair_read_8.simpleName}.bai

  picard CollectHsMetrics I=${pair_read_8} O=${pair_read_8.simpleName}_hs_metrics.txt \
  R=/var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  BAIT_INTERVALS=/var/spool/mail/bait_capture_files/SeqCap_EZ_MedExome/MedExome_hg38_capture_targets.interval \
  TARGET_INTERVALS=/var/spool/mail/bait_capture_files/SeqCap_EZ_MedExome/MedExome_hg38_empirical_targets.interval

  picard CollectInsertSizeMetrics I=${pair_read_8} H=${pair_read_8.simpleName}_histogram.pdf \
  O=${pair_read_8.simpleName}_insert_size.txt M=0.5
  """
}
