
/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output/aligned_sorted/*_processed.bam"
bam_ch = Channel .fromPath( params.bam )

bam_ch.into { bam2_ch; bam3_ch }
tumor_ch = Channel .fromPath (params.tumor )
normal_ch = Channel .fromPath (params.normal )

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


process BaseRecalibrator {
  storeDir "$baseDir/output/aligned_sorted"
  input:
  file pair_read_4 from bam2_ch
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
  file pair_read_6 from bam3_ch
  output:
  file "${pair_read_6.simpleName}.BQSR.bam" into (mutect2_1_ch, mutect2_2_ch)
  script:
  """
  gatk ApplyBQSR \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -I ${pair_read_6} \
  --bqsr-recal-file ${pair_read_5} \
  -O ${pair_read_6.simpleName}.BQSR.bam
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
  file "${x}vs${y}.vcf.gz" into filter_vcf_ch
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

process filter_vcf {
	storeDir "$baseDir/output/mutect2/filtered_vcf"
	input:
	file "${filtered}" from filter_vcf_ch
	output:
	file "${filtered.simpleName}.filtered.vcf" into functotator_ch
	script:
	"""
	gatk FilterMutectCalls \
	-R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
	-V ${filtered}  \
	-O ${filtered.simpleName}.filtered.vcf
	"""
}

process functotator {
  storeDir "$baseDir/output/functotator"
  input:
  file vcf from functotator_ch
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
