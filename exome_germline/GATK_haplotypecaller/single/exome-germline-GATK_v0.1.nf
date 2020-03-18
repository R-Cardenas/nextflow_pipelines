
/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output/aligned_sorted/*_processed.bam"

bam_ch = Channel .fromPath( params.bam )



println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         GATK Germline - single samples
         v0.1
         ===================================



         """
         .stripIndent()

process haplotypeCaller {
  storeDir "$baseDir/output/haplotypeCaller_vcf"
  input:
  file bam from bam_ch
  output:
  file "${bam.simpleName}.vcf.gz" into (haplotype2_ch, index2_ch)
  script:
  """
  gatk HaplotypeCaller \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -I ${bam} \
  -O ${bam.simpleName}.vcf.gz \
  --create-output-variant-index true
  """
}

process IndexFeatureFile {
  storeDir "$baseDir/output/haplotypeCaller_vcf"
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
  storeDir "$baseDir/output/CNNscoreVariants"
  input:
  file vcf from haplotype2_ch
  file idx from idx2_ch
  output:
  file "${vcf.simpleName}.vcf" into ccn2_ch
  script:
  """
  gatk CNNScoreVariants \
  -V ${vcf} \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -O ${vcf.simpleName}.vcf \
  --read-index ${idx}
  """
}

process FilterVariantTranches {
  storeDir "$baseDir/output/FilterVariantTranches"
  input:
  file vcf from ccn2_ch
  output:
  file "${vcf.simpleName}_filtered.vcf" into vcffilter_ch
  script:
  """
  gatk FilterVariantTranches \
  -V ${vcf} \
  --resource /var/spool/mail/hg38/GATK/germline_resource/Homo_sapiens_assembly38.dbsnp138.vcf \
  --resource /var/spool/mail/hg38/GATK/germline_resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --resource /var/spool/mail/hg38/GATK/germline_resource/hapmap_3.3.hg38.vcf.gz \
  --resource /var/spool/mail/hg38/GATK/germline_resource/1000G_omni2.5.hg38.vcf.gz \
  --info-key CNN_1D \
  --snp-tranche 99.95 \
  --indel-tranche 99.4 \
  -O "${vcf.simpleName}_filtered.vcf"
  """
}

process functotator {
  storeDir "$baseDir/output/functotator"
  input:
  file vcf from vcffilter_ch
  output:
  file "${vcf.baseName}.maf" into maf_ch
  script:
  """
  gatk Funcotator \
   -R /var/spool/mail/Homo_sapiens_assembly38.fasta \
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
