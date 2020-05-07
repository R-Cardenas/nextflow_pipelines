
/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output/trim/merge_lanes/*merged.bam"
bam_ch = Channel .fromPath( params.bam )

bam_ch.into { bam2_ch; bam3_ch }

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
  file bam from bam2_ch
  output:
  file "${bam}.table" into table_ch
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
	rm -fr tmp
  """
}

//line 47/55 change
process applyBaseRecalibrator {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort"
  input:
  file bam from bam3_ch
	file "${bam}.table" from table_ch
  output:
  file "${bam.simpleName}.BQSR.bam" into (haplotype_bam_ch, index_ch)
  script:
  """
	mkdir -p tmp
  gatk ApplyBQSR \
  -R $genome_fasta \
  -I ${bam} \
  --bqsr-recal-file ${bam}.table \
  -O ${bam.simpleName}.BQSR.bam \
	--tmp-dir tmp
	rm -fr tmp
  """
}

process bam_index {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort"
  input:
  file pair_read_8 from index_ch
  output:
  file "${pair_read_8}.bai" into haplotype_index_ch

  script:
  """
  samtools index ${pair_read_8}

  """
}
// line 77/86 change
process haplotypeCaller {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/haplotypeCaller"
  input:
	file bam from haplotype_bam_ch
  file "${bam}.bai" from haplotype_index_ch
  output:
  file "${bam.simpleName}.g.vcf.gz" into gatk_combine_ch
  script:
  """
	mkdir -p tmp
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

//change py to bin dir
process combine_gvcf {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/haplotypeCaller"
  input:
  file vcf from gatk_combine_ch.collect()
  output:
  file "${projectname}_combined.g.vcf" into combine_ch
  script:
  """
  python $baseDir/bin/GATK_CombineGVCF.py -V '${vcf}' -O ${projectname}_combined.g.vcf
  """
}

process genotypeVCF {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/haplotypeCaller"
  input:
  file combine from combine_ch
  output:
  file "${combine.simpleName}.genotype.vcf" into (index2_ch, cnn1_ch)
  script:
  """
	mkdir -p tmp
  gatk GenotypeGVCFs \
  -R $genome_fasta \
  -V ${combine} \
  -O ${combine.simpleName}.genotype.vcf \
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
  output:
  file "${vcf.simpleName}.vcf.gz.tbi" into cnn2_ch
  script:
  """
	mkdir -p tmp
  gatk IndexFeatureFile \
    -F ${vcf} \
    -O ${vcf.simpleName}.vcf.gz.tbi \
		--tmp-dir tmp
	rm -fr tmp
  """
}

process CNNscoreVariants {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/haplotypeCaller"
  input:
  file vcf from cnn1_ch
  file "${vcf.simpleName}.vcf.gz.tbi" from cnn2_ch
  output:
  file "${vcf.simpleName}.vcf" into filterVCF_ch
  script:
  """
	mkdir -p tmp
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
	file "${zip}.gz" into merge_ch
	script:
	"""
	bgzip ${zip}
	"""
}

// use bcftools which one?? normal conda version
process merge_vcf {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/GATK_germline_cohort/filtered_vcf"
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
	file "${projectname}_GATK_single_v0.2_filtered.vcf.gz" into index3_ch
	script:
	"""
	cp ${zip2} $baseDir/output/VCF_collect/${projectname}_GATK_cohort.vcf.gz
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
