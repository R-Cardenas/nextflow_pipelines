
/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output/aligned_sorted/*_processed.bam"
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
  storeDir "$baseDir/output/aligned_sorted"
  input:
  file pair_read_4 from bam2_ch
  output:
  file "${pair_read_4.simpleName}_calibration.table" into table_ch
  script:
  """
  gatk BaseRecalibrator \
  -I $pair_read_4 \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  --known-sites /var/spool/mail/hg38/GATK/germline_resource/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites /var/spool/mail/hg38/GATK/germline_resource/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --known-sites /var/spool/mail/hg38/GATK/germline_resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O ${pair_read_4.simpleName}_calibration.table
  """
}

process applyBaseRecalibrator {
  storeDir "$baseDir/output/aligned_sorted"
  input:
  file pair_read_5 from table_ch
  file pair_read_6 from bam3_ch
  output:
  file "${pair_read_6.simpleName}.BQSR.bam" into bam4_ch
  script:
  """
  gatk ApplyBQSR \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -I ${pair_read_6} \
  --bqsr-recal-file ${pair_read_5} \
  -O ${pair_read_6.simpleName}.BQSR.bam
  """
}

process haplotypeCaller {
  storeDir "$baseDir/output/haplotypeCaller"
  input:
  file bam from haplotype_ch
  file index from index_ch
  output:
  file "${bam.simpleName}.g.vcf.gz" into (haplotype2_ch, index2_ch)
  script:
  """
  gatk HaplotypeCaller \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -I ${bam} \
  --read-index ${index} \
  -O ${bam.simpleName}.g.vcf.gz \
  --create-output-variant-index true
  -ERC GVCF
  """
}

process combine_gvcf {
  storeDir "$baseDir/output/haplotypeCaller"
  input:
  file vcf from haplotype2_ch.collectFile()
  output:
  file "project_combined.g.vcf" into combine_ch
  script:
  """
  python GATK_CombineGVCF.py -V '${vcf}' -O project_combined.g.vcf
  """
}

process genotypeVCF {
  storeDir "$baseDir/output/genotypeVCF"
  input:
  file combine from combine_ch
  output:
  file "${combine}.genotype.vcf" into (genotype_ch, haplotype2_ch)
  script:
  """
  gatk GenotypeGVCFs \
  -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa \
  -V ${combing}
  -O ${combine}.genotype.vcf
  """
}

process IndexFeatureFile {
  storeDir "$baseDir/output/haplotypeCaller_vcf"
  input:
  file vcf from genotype_ch
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
