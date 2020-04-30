
/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output/aligned_sorted/*.rename.bam"
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
	-R $genome_fasta \
  --known-sites $GATK_dbsnp138 \
  --known-sites $GATK_1000G \
  --known-sites $GATK_mills \
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
	-R $genome_fasta \
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
  file "${x}.BQSR.bam" from mutect2_1_ch
  file "${y}.BQSR.bam" from mutect2_2_ch
  output:
  val "${x}vs${y}.vcf.gz" into filter_vcf_ch
  script:
  """
  gatk Mutect2 \
	-R $genome_fasta \
  -I ${x}.BQSR.bam \
  -I ${y}.BQSR.bam \
  -normal ${y}_R1_001_BQSR.bam \
  --germline-resource $Mutect2_germline \
  --panel-of-normal $Mutect2_PoN \
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
	-R $genome_fasta \
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
	 -R $genome_fasta \
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
