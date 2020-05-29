// 1. Gather VCF files (done by preceding pipelines)
// 2. Overlap from different tools (eg freebayes / gatk)
// 3. Functotator MAF samtools

//params.freebayes = "$baseDir/output/freebayes/${projectname}_filtered_freebayes.vcf.gz"
params.freebayes = "$baseDir/output/VCF_collect/*.vcf.gz"
params.index = "$baseDir/output/VCF_collect/*.vcf.gz.csi"

vcf_ch = Channel. fromPath (params.freebayes)
vcf_ch.into {vcf2_ch; vcf3_ch}

index_ch = Channel. fromPath (params.index)
index_ch.into {index2_ch; index3_ch}


process vcf_stats {
  storeDir "$baseDir/output/VCF_collect"
  input:
  file vcf from vcf2_ch.collect()
  file index from index2_ch.collect()
  output:
  file "${projectname}_bcf_stats.vchk" into stats_ch
  script:
  """
  bcftools stats -c both ${vcf} > ${projectname}_bcf_stats.vchk
  """
}


process awks_stats {
  storeDir "$baseDir/output/VCF_collect"
  input:
  file vchk from stats_ch
  output:
  file "SN_${projectname}.txt" into sn_ch
  file "ID_${projectname}.txt" into id_ch
  script:
  """
  awk '\$1 == "SN" {print}' ${vchk} > SN_${projectname}.txt
  awk '\$1 == "ID" {print}' ${vchk} > ID_${projectname}.txt
  """
}

//  ## put inside container and chmod +x it the bcftools
process overlap_stats {
  storeDir "$baseDir/output/VCF_collect"
  input:
  file snps from sn_ch
  file id from id_ch
  output:
  file "indel_venn_${projectname}.jpg"
  file "snp_venn_${projectname}.jpg"
  script:
  """
  /bcftools_stat_plot.R \
  -I ${snps} \
  -D ${id} \
  -o snp_venn_${projectname}.jpg \
  -O indel_venn_${projectname}.jpg
  """
}

// add bcftools singularity and targets / regions files for exome
process isec {
  storeDir "$baseDir/output/VCF_collect/plots"
  input:
  file vcf from vcf3_ch.collect()
  file index from index3_ch.collect()
  output:
  file "${projectname}_merged.vcf" into (functotator_ch, somalier_ch)
  script:
  """
  mkdir -p tmp
  bcftools isec -c both ${vcf} -o ${projectname}_merged.vcf -p tmp
  mv tmp/0002.vcf ./${projectname}_merged.vcf
  rm -fr tmp
  """
}
// add to config
process functotator {
  storeDir "$baseDir/output/functotator"
  input:
  file vcf from functotator_ch
  output:
  file "${vcf.baseName}.maf" into maf_ch
  script:
  """
  gatk IndexFeatureFile -F ${vcf}

  gatk Funcotator \
   -R /var/spool/mail/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa \
   -V ${vcf} \
   -O ${vcf.baseName}.maf \
   --output-file-format MAF \
   --data-sources-path /var/spool/mail/GATK_functotator_files/funcotator_dataSources.v1.6.20190124g \
   --ref-version hg38
  """
}

process maf_header {
  storeDir "$baseDir/output/functotator"
  input:
  file maf from maf_ch
  output:
  file "${maf.simpleName}_noheader.maf" into allele_frequency_ch
  script:
  """
  awk '!/\\#/' ${maf} > ${maf.simpleName}_noheader.maf
  """
}

process gnomAD_AF {
  storeDir "$baseDir/output/functotator"
  input:
  file dom from som_ch
  file maf from allele_frequency_ch
  output:
  file "${projectname}_filtered.maf" into maftools_ch
  script:
  """
  chmod +x $baseDir/bin/gnomAD_AF_MAF.R
  $baseDir/bin/gnomAD_AF_MAF.R \
  -f ${maf} \
  -E $AF_group1 \
  -AF $AF1 \
  -x $AF_group2 \
  -y $AF2 \
  -o $projectname \
  -D $DP
  """
}

process maftools {
  storeDir "$baseDir/output/functotator"
  input:
  file maf from maftools_ch
  output:
  file "${projectname}_variant_summary.jpg"
  file "${projectname}_oncoplot.jpg"
  file "${projectname}_titv.jpg"
  file "${projectname}_VAF.jpg"
  script:
  """
  chmod +x $baseDir/bin/maftools.R
  $baseDir/bin/maftools.R \
  -f ${maf}
  -o $projectname
  """
}
