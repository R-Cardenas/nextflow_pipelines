

/*
 * create a channel for fastq pairs
 */
params.bam = "$baseDir/output/trim/merge_lanes/*merged.bam"

Channel
	.fromPath( params.bam )
	.collect ()
	.set {bam_ch}


println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         Freebayes - Germline
         v0.2
         ===================================



         """
         .stripIndent()


process Freebayes {
  storeDir "$baseDir/output/freebayes"
  input:
  file bam from bam_ch
  output:
  file "${projectname}_raw.vcf" into vcf_ch // inset $project name from config file
  script:
  """
  freebayes ${bam} -f $genome_fasta > ${projectname}_raw.vcf
  """
}

process vcf_filter {
  scratch true
  input:
  file vcf from vcf_ch
  output:
  file "${projectname}_filtered_freebayes.vcf" into zip_ch
  script:
  """
  vcffilter -f \
  "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
  ${vcf} > ${projectname}_filtered_freebayes.vcf
  """
}


//needs bcftools
process bam_index {
	storeDir "$baseDir/output/VCF_collect"
	input:
	file vcf from vcf_ch
	output:
	file "${vcf}.csi"
	file "${vcf}.gz"
	script:
	"""
	bgzip ${vcf}
	bcftools index ${vcf}.gz
	"""
}

workflow.onComplete {
	// create log files and record output
	process finish{
		script:
		""" echo ' Pipeline Freebayes germline v0.2 completed
		Project: $projectname
		Time: ${nextflow.timestamp}

	 Freebayes - completed
	 (Reference - $genome_fasta)
	 VCFTools - completed
	 (Parameters: "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1")
	 Bgzip - completed
	 BcfTools Index - completed
	 ' >> $baseDir/${projectname}_log.txt

	 mail -s "Freebayes successful" < $baseDir/${projectname}_log.txt
	 """
	}
}

workflow.onError {
	process finish_error{
		script:
		""" echo 'Pipeline Freebayes FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}

		Error:
		${workflow.errorMessage}' {input} >> $baseDir/${projectname}_error.txt

	  mail -s "cgpMAP successful" {input} < $baseDir/${projectname}_error.txt
	  """
	}
}
