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

// name in nextflow config
process vcf_stats {
  storeDir "$baseDir/output/VCF_collect"
  input:
  file vcf from vcf2_ch.collect()
  file index from index2_ch.collect()
  output:
  file "${projectname}_bcf_stats.vchk" into stats_ch
  script:
  """
  ls ${index}
  bcftools stats -c none ${vcf} > ${projectname}_bcf_stats.vchk
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
  awk '$1 == "SN" {print}' ${vchk} > SN_${projectname}.txt
  awk '$1 == "ID" {print}' ${vchk} > ID_${projectname}.txt
  """
}

process plots {
  storeDir "$baseDir/output/VCF_collect"
}


// add bcftools singularity and targets / regions files for exome
process isec {
  storeDir "$baseDir/output/VCF_collect/plots"
  input:
  file vcf from vcf3_ch.collect()
  file index from index3_ch.collect()
  output:
  file "${projectname}_merged.vcf.gz" into functotator_ch
  script:
  """
  ls ${index}
  bcftools isec -c none ${vcf} -O z -o ${projectname}_merged.vcf.gz -p tmp
  mv tmp/0002.vcf.gz ./${projectname}_merged.vcf.gz
  rm -fr tmp
  """
}

// add to config
process functotator {
  storeDir "$baseDir/output/hg38_decoy/GATK_germline_single/gatk_germline_single/functotator"
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
