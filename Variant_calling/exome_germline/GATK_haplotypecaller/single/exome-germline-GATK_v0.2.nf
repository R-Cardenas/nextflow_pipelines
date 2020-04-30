
/*
 * create a channel for bam files produced by Pipeline GATK germline (single)_processing pipeline
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
  storeDir "$baseDir/output/GATK_germline_single"
  input:
  file bam from bam2_ch
  output:
  file "${bam.simpleName}_calibration.table" into table_ch
  script:
  """
  gatk BaseRecalibrator \
  -I ${bam} \
  -R $genome_fasta \
  --known-sites $GATK_dbsnp138 \
  --known-sites $GATK_1000G \
  --known-sites $GATK_mills \
  -O ${bam.simpleName}_calibration.table
  """
}

process applyBaseRecalibrator {
  storeDir "$baseDir/output/GATK_germline_single"
  input:
  file table from table_ch
  file bam from bam3_ch
  output:
  file "${bam.simpleName}.BQSR.bam" into (haplotype_ch, index_ch)
  script:
  """
  gatk ApplyBQSR \
  -R $genome_fasta \
  -I ${bam} \
  --bqsr-recal-file ${table} \
  -O ${bam.simpleName}.BQSR.bam
  """
  }

process bam_index {
  storeDir "$baseDir/output/GATK_germline_single"
  input:
  file index from index_ch
  output:
  file "${index}.bai" into index2_ch

  script:
  """
  samtools index ${index}

  """
}

process haplotypeCaller {
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller"
  input:
  file bam from haplotype_ch
  file bai from index2_ch
  output:
  file "${bam.simpleName}.vcf.gz" into (haplotype2_ch, index3_ch)
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
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller"
  input:
  file vcf from index3_ch
  output:
  file "${vcf.simpleName}.vcf.gz.tbi" into cnn_ch
  script:
  """
  gatk IndexFeatureFile \
    -F ${vcf} \
    -O ${vcf.simpleName}.vcf.gz.tbi
  """
}

process CNNscoreVariants {
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller"
  input:
  file vcf from haplotype2_ch
  file idx from cnn_ch
  output:
  file "${vcf.simpleName}.vcf" into filterVCF_ch
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
  storeDir "$baseDir/output/GATK_germline_cohort/filtered_vcf"
	input:
	file vcf from filterVCF_ch
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
  storeDir "$baseDir/output/GATK_germline_cohort/filtered_vcf"
	input:
	file zip from zip_ch
	output:
	file "${zip}.gz" into (merge_ch, index4_ch)
	script:
	"""
	bgzip ${zip}
	"""
}

//needs bcftools
process bcf_index1 {
  storeDir "$baseDir/output/GATK_germline_cohort/filtered_vcf"
	input:
	file vcf from index4_ch
	output:
	file "${vcf}.csi" into merge2_ch
	script:
	"""
	bcftools index ${vcf}
	"""
}

// use bcftools which one?? normal conda version
process merge_vcf {
  storeDir "$baseDir/output/GATK_germline_cohort/filtered_vcf"
	input:
	file vcf2 from merge_ch.collect()
	file index from merge2_ch.collect()
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
process bcf_index {
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

workflow.onComplete {
	// create log files and record output
	process finish{
		script:
		""" echo ' Pipeline GATK germline (single) v0.2 completed
		Project: $projectname
		Time: ${nextflow.timestamp}
		BaseRecalibrator - completed
		(Reference: $genome_fasta
		 known sites: $GATK_dbsnp138
		 known sites: $GATK_1000G
		 known sites: $GATK_mills)

		applyBaseRecalibrator - completed
		bam index  - completed
		haplotypeCaller - completed

	 ' >> $baseDir/${projectname}_log.txt

	 mail -s "GATK germline (single) successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_log.txt
	 """
	}
}

workflow.onError {
	process finish_error{
		script:
		"""
		echo 'Pipeline GATK germline (single) v0.2 FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}

		Error:
		${workflow.errorMessage}' >> $baseDir/${projectname}_error.txt

	  mail -s "Pipeline GATK germline (single) unsuccessful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_error.txt
	  """
	}
}
