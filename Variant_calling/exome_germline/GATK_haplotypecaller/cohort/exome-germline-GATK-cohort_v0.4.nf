
/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output/trim/merge_lanes/*merged.bam"
bam_ch = Channel .fromPath( params.bam )


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

process BaseRecalibrator {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort"
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
  storeDir "$baseDir/output/GATK_germline_cohort/haplotypeCaller"
  input:
	file bam from haplotype_bam_ch
  output:
  file "${bam.simpleName}.g.vcf.gz" into (vcf_index, gatk_combine_ch)
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
  -ERC GVCF \
	--tmp-dir tmp
	rm -fr tmp
  """
}

process vcf_index{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/haplotypeCaller"
	input:
  file vcf from vcf_index
  output:
	file "${vcf}.*" into vcf_index2_ch
	script:
	"""
	mkdir -p tmp
  gatk IndexFeatureFile \
    -F ${vcf} \
		--tmp-dir tmp
	"""
}

//change py to bin dir
process combine_gvcf {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/haplotypeCaller"
  input:
  file vcf from gatk_combine_ch.collect()
	file index from vcf_index2_ch.collect()
  output:
  file "${projectname}_combined.g.vcf.gz" into combine_ch
  script:
  """
	mkdir -p tmp
  python $baseDir/bin/GATK_CombineGVCF.py -V '${vcf}' -O ${projectname}_combined.g.vcf.gz -R $genome_fasta
  """
}

process genotypeVCF {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/haplotypeCaller"
  input:
  file combine from combine_ch
  output:
  file "${combine.simpleName}.genotype.vcf.gz" into index2_ch
	file "${vcf.simpleName}.vcf.gz.tbi" into index10_ch
  script:
  """
	mkdir -p tmp
  gatk IndexFeatureFile \
    -F ${combine} \
    -O ${combine.simpleName}.vcf.gz.tbi \
		--tmp-dir tmp

  gatk GenotypeGVCFs \
  -R $genome_fasta \
  -V ${combine} \
  -O ${combine.simpleName}.genotype.vcf.gz \
	--tmp-dir tmp
	rm -fr tmp
  """
}

process IndexFeatureFile {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/haplotypeCaller"
  input:
  file vcf from index2_ch
	file index from index10_ch
  output:
	file "${vcf.simpleName}.vcf" into filterVCF_ch
  script:
  """
	mkdir -p tmp
  gatk IndexFeatureFile \
    -F ${vcf} \
    -O ${vcf.simpleName}.vcf.gz.tbi \
		--tmp-dir tmp

	gatk CNNScoreVariants \
	  -V ${vcf} \
	  -R $genome_fasta \
	  -O ${vcf.simpleName}.vcf \
	  --read-index ${vcf.simpleName}.vcf.gz.tbi \
		--tmp-dir tmp
	rm -fr tmp
  """
}

process FilterVariantTranches {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/filtered_vcf"
  input:
  file vcf from filterVCF_ch
  output:
  file "${vcf.simpleName}_filtered.vcf" into zip_ch
  script:
  """
	mkdir -p tmp
  gatk FilterVariantTranches \
  -V ${vcf} \
	--resource $GATK_dbsnp138 \
  --resource $GATK_1000G \
  --resource $GATK_mills \
	--resource $GATK_hapmap \
  --info-key CNN_1D \
  --snp-tranche 99.95 \
  --indel-tranche 99.4 \
  -O "${vcf.simpleName}_filtered.vcf" \
	--tmp-dir tmp
	rm -fr tmp
  """
}

// needs samttools
process zip {
  storeDir "$baseDir/output/GATK_germline_cohort/filtered_vcf"
	input:
	file zip from zip_ch
	output:
	file "${zip}.gz" into collect_ch
	script:
	"""
	bgzip ${zip}
	"""
}

process collect_vcf {
	storeDir "$baseDir/output/VCF_collect"
	input:
	file zip2 from collect_ch
	output:
	file "${projectname}_GATK_single_v0.2_filtered.vcf.gz" into index3_ch
	script:
	"""
	mkdir -p $baseDir/output/VCF_collect
	cp ${zip2} $baseDir/output/VCF_collect/${projectname}_GATK_single_v0.2_filtered.vcf.gz
	"""
}

//needs bcftools
process bcf_index {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/VCF_collect"
	input:
	file vcf from index3_ch
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
		""" echo ' Pipeline GATK germline (cohort) v0.3 completed
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

	 mail -s "GATK germline (cohort) successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_log.txt
	 """
	}
}

workflow.onError {
	process finish_error{
		script:
		"""
		echo 'Pipeline cgpmap v0.3 FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}

		Error:
		${workflow.errorMessage}' >> $baseDir/${projectname}_error.txt

	  mail -s "cgpMAP successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_error.txt
	  """
	}
}
