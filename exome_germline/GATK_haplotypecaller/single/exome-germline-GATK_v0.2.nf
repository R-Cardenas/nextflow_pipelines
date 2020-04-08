
/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output/aligned_sorted/*.rename.bam"
bam_ch = Channel .fromPath( params.bam )
bam_ch.into { bam2_ch; bam3_ch }

println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         GATK Germline - HaplotypeCaller
                            Single samples
         v0.2
         ===================================



         """
         .stripIndent()

process BaseRecalibrator {
  storeDir "$baseDir/output/GATK_germline_single/aligned_sorted"
  input:
  file pair_read_4 from bam2_ch
  output:
  file "${pair_read_4.simpleName}_calibration.table" into table_ch
  script:
  """
  gatk BaseRecalibrator \
  -I ${pair_read_4} \
  -R $genome_fasta \
  --known-sites $GATK_dbsnp138 \
  --known-sites $GATK_1000G \
  --known-sites $GATK_mills \
  -O ${pair_read_4.simpleName}_calibration.table
  """
}

process applyBaseRecalibrator {
  storeDir "$baseDir/output/GATK_germline_single/aligned_sorted"
  input:
  file pair_read_5 from table_ch
  file pair_read_6 from bam3_ch
  output:
  file "${pair_read_6.simpleName}.BQSR.bam" into (bam4_ch, index_66ch)
  script:
  """
  gatk ApplyBQSR \
  -R $genome_fasta \
  -I ${pair_read_6} \
  --bqsr-recal-file ${pair_read_5} \
  -O ${pair_read_6.simpleName}.BQSR.bam
  """
  }

process bam_index {
  storeDir "$baseDir/output/GATK_germline_single/aligned_sorted"
  input:
  file pair_read_8 from index_66ch
  output:
  file "${pair_read_8}.bai" into index100_ch

  script:
  """
  samtools index ${pair_read_8}

  """
}

process haplotypeCaller {
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller_vcf"
  input:
  file bam from bam4_ch
  file bai from index100_ch
  output:
  file "${bam.simpleName}.vcf.gz" into (haplotype2_ch, index2_ch)
  script:
  """
  gatk HaplotypeCaller \
  -R $genome_fasta \
  -I ${bam} \
  -O ${bam.simpleName}.vcf.gz \
  --create-output-variant-index true
  """
}

process IndexFeatureFile {
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller_vcf"
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
  scratch true
  input:
  file vcf from haplotype2_ch
  file idx from idx2_ch
  output:
  file "${vcf.simpleName}.vcf" into ccn2_ch
  script:
  """
  gatk CNNScoreVariants \
  -V ${vcf} \
  -R $genome_fasta \
	-O ${vcf.simpleName}.vcf \
	--read-index ${idx}
	"""
}

process FilterVariantTranches {
	scratch true
	input:
	file vcf from ccn2_ch
	output:
	file "${vcf.simpleName}_filtered.vcf" into zip_ch
	script:
	"""
	gatk FilterVariantTranches \
	-V ${vcf} \
	--resource $GATK_dbsnp138 \
  --resource $GATK_1000G \
  --resource $GATK_mills \
	--resource $GATK_hapmap \
	--info-key CNN_1D \
	--snp-tranche 99.95 \
	--indel-tranche 99.4 \
	-O "${vcf.simpleName}_filtered.vcf"
	"""
}

// needs samttools
process zip {
  storeDir "$baseDir/output/GATK_germline_single/FilteredVCF"
	input:
	file zip from zip_ch
	output:
	file "${zip}.gz" into merge_ch
	script:
	"""
	bgzip ${zip}
	"""
}

// use bcftools which one?? normal conda version
process merge_vcf {
  storeDir "$baseDir/output/GATK_germline_single/FilteredVCF"
	input:
	file vcf2 from merge_ch.collect()
	output:
	file "${projectname}_GATK_single_v0.2_filtered.vcf.gz" into collect_ch
	script:
	"""
	bcftools merge -m all ${vcf2} -O z -o ${projectname}_GATK_single_v0.2_filtered.vcf.gz
	"""
}

process collect_vcf {
	storeDir "$baseDir/output/VCF_collect"
	input:
	file zip2 from collect_ch
	output:
	file "${projectname}_GATK_single_v0.2_filtered.vcf.gz" into index101_ch
	script:
	"""
	cp ${zip2} $baseDir/output/VCF_collect/${projectname}_GATK_single_v0.2_filtered.vcf.gz
	"""
}

//needs bcftools
process index {
	storeDir "$baseDir/output/VCF_collect"
	input:
	file vcf from index101_ch
	output:
	file "${vcf}.csi"
	script:
	"""
	bcftools index ${vcf}
	"""
}
