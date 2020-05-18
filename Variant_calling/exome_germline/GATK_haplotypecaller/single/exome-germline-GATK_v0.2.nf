
/*
 * create a channel for bam files produced by Pipeline GATK germline (single)_processing pipeline
 */
params.bam = "$baseDir/output/BAM/merged_lanes/*.rmd.bam"
bam_ch = Channel .fromPath( params.bam )


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
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller"
  input:
  file bam from bam_ch
  output:
  file "${bam.simpleName}.BQSR.bam" into haplotype_bam_ch
	file "${bam}.table"
  script:
  """
	mkdir -p tmp
  gatk BaseRecalibrator \
	-I ${bam} \
  -R $genome_fasta \
  --known-sites $GATK_dbsnp138 \
  --known-sites $GATK_1000G \
  --known-sites $GATK_mills \
  -O ${bam}.table \
	--tmp-dir tmp

	gatk ApplyBQSR \
  -R $genome_fasta \
  -I ${bam} \
  --bqsr-recal-file ${bam}.table \
  -O ${bam.simpleName}.BQSR.bam \
	--tmp-dir tmp
	rm -fr tmp

  """
}

// line 77/86 change
process haplotypeCaller {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller"
  input:
	file bam from haplotype_bam_ch
  output:
  file "${bam.simpleName}.g.vcf.gz" into haplotype2_ch
	file "${bam}.bai"
  script:
  """
	mkdir -p tmp
	gatk BuildBamIndex \
	-I ${bam} \
	-O ${bam}.bai \
	--TMP_DIR tmp

  gatk HaplotypeCaller \
  -R $genome_fasta \
  -I ${bam} \
	--read-index ${bam}.bai \
  -O ${bam.simpleName}.g.vcf.gz \
  --create-output-variant-index true \
	--intervals $haplotypecaller_bed \
  -ERC NONE \
	--tmp-dir tmp
	rm -fr tmp
  """
}


process CNNscoreVariants {
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller"
  input:
  file vcf from haplotype2_ch
  output:
  file "${vcf.simpleName}_CNN.g.vcf.gz" into filterVCF_ch
  script:
  """
	mkdir -p tmp
	gatk IndexFeatureFile \
		-F ${vcf} \
		--tmp-dir tmp

  gatk CNNScoreVariants \
  -V ${vcf} \
  -R $genome_fasta \
	-O ${vcf.simpleName}_CNN.g.vcf.gz \
	--tmp-dir tmp
	"""
}

process FilterVariantTranches {
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller"
	input:
	file vcf from filterVCF_ch
	output:
	file "${vcf.simpleName}_filtered.g.vcf" into zip_ch
	script:
	"""
	mkdir -p tmp
	gatk IndexFeatureFile \
		-F ${vcf} \
		--tmp-dir tmp

	gatk FilterVariantTranches \
	-V ${vcf} \
	--info-key CNN_1D \
	--resource $GATK_dbsnp138 \
  --resource $GATK_1000G \
  --resource $GATK_mills \
	--resource $GATK_hapmap \
	--snp-tranche 99.95 \
	--indel-tranche 99.4 \
	-O ${vcf.simpleName}_filtered.g.vcf \
	--tmp-dir tmp
	"""
}

// needs samttools
process zip {
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller"
	input:
	file zip from zip_ch
	output:
	file "${zip}.gz" into merge2_ch
	file "${zip}.gz.csi" into merge3_ch
	script:
	"""
	bgzip ${zip}
	bcftools index ${zip}.gz
	"""
}

// use bcftools which one?? normal conda version
process merge_vcf {
  storeDir "$baseDir/output/GATK_germline_single/haplotypeCaller"
	input:
	file vcf2 from merge2_ch.collect()
	file index2 from merge3_ch.collect()
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
