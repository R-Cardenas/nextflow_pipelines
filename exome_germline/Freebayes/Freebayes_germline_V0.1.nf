/*
 * create a channel for fastq pairs
 */
params.bam = "$baseDir/output/aligned_sorted/*.rename.bam"

bam_ch = Channel .fromPath( params.bam )


println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         Freebayes - Germline
         v0.1
         ===================================



         """
         .stripIndent()

process Freebayes {
  storeDir "$baseDir/output/freebayes"
  input:
  file bam from bam_ch.collectFile()
  output:
  file "projectname_raw.vcf" into vcf_ch // inset $project name from config file
  script:
  """
  freebayes ${bam} -f /var/spool/mail/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa > projectname_raw.vcf
  """
}

process vcf_filter {
  storeDir "$baseDir/output/freebayes"
  input:
  file vcf from vcf_ch
  output:
  file "projectname_processed.vcf"
  script:
  """
  vcffilter -f \
  "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
  ${vcf} > projectname_processed.vcf
  """
}
