
/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output/aligned_sorted/*.rename.bam"
bam_ch = Channel .fromPath( params.bam )

bam_ch.into { bam2_ch; bam3_ch }
tumor_ch = Channel .from (params.tumor )
normal_ch = Channel .from (params.normal )

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
  storeDir "$baseDir/output/mutect2/aligned_sorted"
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
  storeDir "$baseDir/output/mutect2/aligned_sorted"
  input:
  file pair_read_5 from table_ch
  file pair_read_6 from bam3_ch
  output:
  file "${pair_read_6.simpleName}.BQSR.bam" into (mutect2_1_ch, mutect2_2_ch, pileup_ch)
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
  storeDir "$baseDir/output/mutect2/mutect2"
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
	-R $genome_fasta \
  -I ${x}_R1_001_BQSR.bam \
  -I ${y}_R1_001_BQSR.bam \
  -normal ${y}_R1_001_BQSR.bam \
  --germline-resource $Mutect2_germline \
  --panel-of-normals $Mutect2_PoN \
  -O ${x}vs${y}.vcf.gz
    """
}

process pileup_summary{
  storeDir "$baseDir/output/mutect2/mutect2"
	input:
	file bam from pileup_ch
	output:
	file "${bam.simpleName}.getpileupsummaries.table" into contamination_ch
	script:
	"""
	gatk GetPileupSummaries \
	-I ${bam} \
	-V $Mutect2_germline \
	-O ${bam.simpleName}.getpileupsummaries.table
	"""
}

process calculate_contamination{
	storeDir "$baseDir/output/mutect2/mutect2"
	input:
	file table from contamination_ch
	output:
	file "${table.simpleName}_calculatecontamination.table" into filter_vcf2_ch
	script:
	"""
	gatk CalculateContamination \
	-I ${table} \
	-O ${table.simpleName}_calculatecontamination.table
	"""
}


process filter_vcf {
	storeDir "$baseDir/output/mutect2/mutect2/filtered_vcf"
	input:
	file vcf from filter_vcf_ch
	file table from filter_vcf2_ch
	output:
	file "${vcf.simpleName}.filtered.vcf" into zip_ch
	script:
	"""
	gatk FilterMutectCalls \
	-R $genome_fasta \
	-V ${vcf}  \
	--contamination-table ${table}
	-O ${vcf.simpleName}.filtered.vcf
	"""
}

// needs samttools
process zip {
  storeDir "$baseDir/output/mutect2/GATK_germline_cohort/FilteredVCF"
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
  storeDir "$baseDir/output/mutect2/GATK_germline_cohort/FilteredVCF"
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
	storeDir "$baseDir/output/mutect2/VCF_collect"
	input:
	file zip2 from collect_ch
	output:
	file "${projectname}_GATK_single_v0.2_filtered.vcf.gz" into index3_ch
	script:
	"""
	cp ${zip2} $baseDir/output/mutect2/VCF_collect/${projectname}_GATK_single_v0.2_filtered.vcf.gz
	"""
}

//needs bcftools
process bcf_index {
	storeDir "$baseDir/output/mutect2/VCF_collect"
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
		Mutect2- completed
		(germline-resource: $Mutect2_germline \
		 panel-of-normal: $Mutect2_PoN )

		Pile up summary - completed
		Contamination calculation - completed
		Filter Mutect VCF - completed


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
