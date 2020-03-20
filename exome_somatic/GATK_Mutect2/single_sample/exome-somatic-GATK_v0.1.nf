
/*
 * create a channel for fastq pairs
 */
params.read1 = "$baseDir/input/*_R1_001.fastq.gz"
params.read2 = "$baseDir/input/*_R2_001.fastq.gz"

read1_ch = Channel .fromPath( params.read1 )
read2_ch = Channel .fromPath( params.read2 )

tumor_ch = Channel .from ( params.tumor)
normal_ch = Channel .from ( params.normal )



println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         GATK Somatic: Mutect2
				 Single samples
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
  file "${read1.simpleName}.fp.fq.gz" into fastp1_ch
  file "${read2.simpleName}.fp.fq.gz" into fastp2_ch
  file "${read1.simpleName}_fastp.json"
  file "${read1.simpleName}_fastp.html"
  script:
  """
  fastp -i ${read1} -I ${read2} -o ${read1.simpleName}.fp.fq.gz -O ${read2.simpleName}.fp.fq.gz \
  --json ${read1.simpleName}_fastp.json \
  --html ${read1.simpleName}_fastp.html
  """
}

process cgpMAP {
  storeDir "$baseDir/output/cgpMAP"
  input:
  file read11 from fastp1_ch
  file read12 from fastp2_ch

  output:
  file "${read11.simpleName}.bam" into cgp_ch

  script:
  """
  ds-cgpmap.pl  \
  -r /var/spool/mail/hg38/UCSC/WholeGenome_hg38.tar.gz \
  -i /var/spool/mail/hg38/UCSC/BWAIndex_hg38.tar.gz \
  -s ${read11.simpleName} \
  -t 6 \
  -outdir $baseDir/output/cgpMAP \
  ${read11} ${read12}
  """
}

process sam_sort {
  storeDir "$baseDir"
  input:
  file pair_read_2 from cgp_ch
  output:
  file "${pair_read_2.baseName}" into sort_ch
  script:
  """
  samtools sort $pair_read_2 -o ${pair_read_2.baseName}
  """
}


process picard_pcr_removal {
  storeDir "$baseDir/output/picard_pcr_removal"
  input:
  file pair_read_3 from sort_ch
  output:
  file "${pair_read_3.baseName}.bam" into rename_ch
  file "${pair_read_3.baseName}_hs_metrics.txt"
  file "${pair_read_3.baseName}_histogram.pdf"
  file "${pair_read_3.baseName}_insert_size.txt"
  script:
  """
  picard MarkDuplicates I=${pair_read_3} O=${pair_read_3.baseName}.bam M=${pair_read_3.baseName}.log

  picard CollectHsMetrics I=${pair_read_3} O=${pair_read_3.baseName}_hs_metrics.txt \
  R=/var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  BAIT_INTERVALS=/var/spool/mail/bait_capture_files/SeqCap_EZ_MedExome/MedExome_hg38_capture_targets.interval \
  TARGET_INTERVALS=/var/spool/mail/bait_capture_files/SeqCap_EZ_MedExome/MedExome_hg38_empirical_targets.interval

  picard CollectInsertSizeMetrics I=${pair_read_3} H=${pair_read_3.baseName}_histogram.pdf \
  O=${pair_read_3.baseName}_insert_size.txt M=0.5

  """
}

process rename_bam {
  storeDir "$baseDir"
  input:
  file pair_read_7 from rename_ch
  output:
  file "${pair_read_7.baseName}_f.bam" into (recal_1_ch, recal_2_ch)
  script:
  """
  picard AddOrReplaceReadGroups I=$pair_read_7 O=${pair_read_7.baseName}_f.bam \
  RGID=rename RGLB=${pair_read_7.baseName} RGPL=illumina RGPU=unit1 RGSM=${pair_read_7.baseName}
  """
}

process BaseRecalibrator {
  storeDir "$baseDir/output/BaseRecalibrator"
  input:
  file pair_read_4 from recal_1_ch
  output:
  file "${pair_read_4.baseName}_calibration.table" into table_ch
  script:
  """
  gatk BaseRecalibrator \
  -I $pair_read_4 \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  --known-sites /var/spool/mail/hg38/GATK/germline_resource/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites /var/spool/mail/hg38/GATK/germline_resource/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --known-sites /var/spool/mail/hg38/GATK/germline_resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O ${pair_read_4.baseName}_calibration.table
  """
}

process applyBaseRecalibrator {
  storeDir "$baseDir/output/BaseRecalibrator"
  input:
  file pair_read_5 from table_ch
  file pair_read_6 from recal_2_ch
  output:
  file "${pair_read_6.baseName}_BQSR.bam" into (bam_ch, mutect2_1_ch, mutect2_2_ch)
  script:
  """
  gatk ApplyBQSR \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -I ${pair_read_6} \
  --bqsr-recal-file ${pair_read_5} \
  -O ${pair_read_6.baseName}_BQSR.bam
  """
}


process mutect2 {
  storeDir "$baseDir/output/mutect2"
  input:
  val x from tumor_ch
  val y from normal_ch
  file "${x}_R1_001_BQSR.bam" from mutect2_1_ch
  file "${y}_R1_001_BQSR.bam" from mutect2_2_ch
  output:
  file "${x}vs${y}.vcf.gz"
  script:
  """
  gatk Mutect2 \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -I ${x}_R1_001_BQSR.bam \
  -I ${y}_R1_001_BQSR.bam \
  -normal ${y}_R1_001_BQSR.bam \
  --germline-resource /var/spool/mail/hg38/GATK/GATK_pon_germline_resource_hg38/somatic-hg38_af-only-gnomad.hg38.vcf \
  --panel-of-normal /var/spool/mail/hg38/GATK/GATK_pon_germline_resource_hg38/somatic-hg38_1000g_pon.hg38.vcf \
  -O ${x}vs${y}.vcf.gz
    """
}
