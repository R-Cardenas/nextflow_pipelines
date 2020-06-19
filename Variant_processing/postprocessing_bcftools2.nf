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
  bcftools stats -c indels ${vcf} > ${projectname}_bcf_stats.vchk
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
  file "${projectname}_merged.vcf" into split_ch
  script:
  """
  mkdir -p tmp
  bcftools isec -c indels ${vcf} -o ${projectname}_merged.vcf -p tmp
  mv tmp/0002.vcf ./${projectname}_merged.vcf
  rm -fr tmp
  """
}

process split_vcf {
  storeDir "$baseDir/output/VCF_collect/split_vcf"
  input:
  file vcf from split_ch
  output:
  file "*.vcf.gz" into (vep_ch, maf_ch)
  script:
  """
  for file in *.vcf*; do
    for sample in `bcftools query -l \$file`; do
      bcftools view -c1 -Oz -s \$sample -o \${file/.vcf*/.\$sample.vcf.gz} \$file
    done
  done
  """
}

// you may want to repeat this  but to create the VCF files also
process VEP {
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP"
  input:
  file vcf from vep_ch.flatten()
  output:
  file "${vcf.simpleName}_VEP.txt"
  script:
  """
  /ensembl-vep/vep -i ${vcf} \
  --dir /var/spool/mail/VEP/.vep \
  -o ${vcf.simpleName}_VEP.txt \
  --cache homo_sapiens \
  --force_overwrite \
  --sift b \
  --polyphen b \
  --variant_class \
  --regulatory \
  --af_gnomad \
  --domains \
  --tab \
  --show_ref_allele \
  --no_headers \
  --verbose
  """
}

process functotator {
  storeDir "$baseDir/output/VCF_collect/split_vcf/functotator"
  input:
  file vcf from maf_ch.flatten()
  output:
  file "${vcf.baseName}.maf" into header_ch
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
  storeDir "$baseDir/output/split_vcf/functotator"
  input:
  file maf from header_ch
  output:
  file "${maf.simpleName}_noheader.maf" into allele_frequency_ch
  script:
  """
  awk '!/\\#/' ${maf} > ${maf.simpleName}_noheader.maf
  """
}
