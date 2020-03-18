
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
         Germline
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
  storeDir "$baseDir/output/aligned_sorted"
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
  storeDir "$baseDir/output/aligned_sorted"
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
  file "${pair_read_7.simpleName}.rename.bam" into (recal_1_ch, recal_2_ch)
  script:
  """
  picard AddOrReplaceReadGroups I=$pair_read_7 O=${pair_read_7.simpleName}.rename.bam \
  RGID=rename RGLB=${pair_read_7.simpleName} RGPL=illumina RGPU=unit1 RGSM=${pair_read_7.simpleName}
  """
}

process BaseRecalibrator {
  storeDir "$baseDir/output/aligned_sorted"
  input:
  file pair_read_4 from recal_1_ch
  output:
  file "${pair_read_4.simpleName}_calibration.table" into table_ch
  script:
  """
  gatk BaseRecalibrator \
  -I $pair_read_4 \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  --known-sites /var/spool/mail/hg38/GATK/germline_resource/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites /var/spool/mail/hg38/GATK/germline_resource/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --known-sites /var/spool/mail/hg38/GATK/germline_resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O ${pair_read_4.simpleName}_calibration.table
  """
}

process applyBaseRecalibrator {
  storeDir "$baseDir/output/aligned_sorted"
  input:
  file pair_read_5 from table_ch
  file pair_read_6 from recal_2_ch
  output:
  file "${pair_read_6.simpleName}.BQSR.bam" into (bam_ch, haplotype_ch)
  script:
  """
  gatk ApplyBQSR \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -I ${pair_read_6} \
  --bqsr-recal-file ${pair_read_5} \
  -O ${pair_read_6.simpleName}.BQSR.bam
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


process haplotypeCaller {
  storeDir "$baseDir/output/haplotypeCaller_vcf"
  input:
  file bam from haplotype_ch
  file index from index_ch
  output:
  file "${bam.simpleName}.vcf.gz" into (haplotype2_ch, index2_ch)
  script:
  """
  gatk HaplotypeCaller \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -I ${bam} \
  --read-index ${index} \
  -O ${bam.simpleName}.vcf.gz \
  --create-output-variant-index true
  """
}

process IndexFeatureFile {
  storeDir "$baseDir/output/haplotypeCaller_vcf"
  input:
  file vcf from index2_ch
  output:
  file "${vcf.simpleName}.vcf.gz.tbi" into idx2_ch
  script:
  """
  gatk IndexFeatureFile \
    -F ${vcf} \
    -O ${vcf.simpleName}.vcf.gz.tbi
  """
}

process CNNscoreVariants {
  storeDir "$baseDir/output/CNNscoreVariants"
  input:
  file vcf from haplotype2_ch
  file idx from idx2_ch
  output:
  file "${vcf.simpleName}.vcf" into ccn2_ch
  script:
  """
  gatk CNNScoreVariants \
  -V ${vcf} \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -O ${vcf.simpleName}.vcf \
  --read-index ${idx}
  """
}

process FilterVariantTranches {
  storeDir "$baseDir/output/FilterVariantTranches"
  input:
  file vcf from ccn2_ch
  output:
  file "${vcf.simpleName}_filtered.vcf" into vcffilter_ch
  script:
  """
  gatk FilterVariantTranches \
  -V ${vcf} \
  --resource /var/spool/mail/hg38/GATK/germline_resource/Homo_sapiens_assembly38.dbsnp138.vcf \
  --resource /var/spool/mail/hg38/GATK/germline_resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --resource /var/spool/mail/hg38/GATK/germline_resource/hapmap_3.3.hg38.vcf.gz \
  --resource /var/spool/mail/hg38/GATK/germline_resource/1000G_omni2.5.hg38.vcf.gz \
  --info-key CNN_1D \
  --snp-tranche 99.95 \
  --indel-tranche 99.4 \
  -O "${vcf.simpleName}_filtered.vcf"
  """
}

process functotator {
  storeDir "$baseDir/output/functotator"
  input:
  file vcf from vcffilter_ch
  output:
  file "${vcf.baseName}.maf" into maf_ch
  script:
  """
  gatk Funcotator \
   -R /var/spool/mail/Homo_sapiens_assembly38.fasta \
   -V ${vcf} \
   -O ${vcf.baseName}.maf \
   --output-file-format MAF \
   --data-sources-path /var/spool/mail/GATK_functotator_files/funcotator_dataSources.v1.6.20190124g \
   --ref-version hg38
  """
}


process process_maf {
  storeDir "$baseDir/output/functotator/processed"
  input:
  file maf from maf_ch
  output:
  file "${maf.baseName}_nohead.maf" into maf2_ch
  script:
  """
  awk '!/#/ {print}' ${maf} > ${maf.baseName}_nohead.maf
  """
}

process multiqc{
  storeDir "$baseDir/output/multiQC"
  input:
  val vcf from maf2_ch.collectFile()
  output:
  file "multiqc_report.html"

  script:
  """
  for i in ${vcf}
  do
    cat \$i >> processed_samples_list.txt
  done
  multiqc $baseDir
  """
}
