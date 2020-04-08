/*
 * create a channel for fastq pairs
 */
params.bam = "$baseDir/output/aligned_sorted/*.rename.bam"

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

// needs samtools
process zip {
  storeDir "$baseDir/output/freebayes"
	input:
	file zip from zip_ch
	output:
	file "${projectname}_filtered_freebaye.vcf.gz" into collect_ch
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
	file "${projectname}_filtered_freebayes.vcf.gz" into index101_ch
	script:
	"""
	cp ${zip2} $baseDir/output/VCF_collect/${projectname}_filtered_freebayes.vcf.gz
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
