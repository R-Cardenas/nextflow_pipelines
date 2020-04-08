// 1. Gather VCF files (done by preceding pipelines)
// 2. Overlap from different tools (eg freebayes / gatk)
// 3. Functotator MAF samtools

//params.freebayes = "$baseDir/output/freebayes/${projectname}_filtered_freebayes.vcf.gz"
params.freebayes = "$baseDir/output/VCF_collect/*.vcf.gz"

vcf_ch = Channel. fromPath (params.freebayes)

vcf_ch.into {vcf2_ch; vcf3_ch}

// name in nextflow config
process vcf_stats {
  storeDir: "$baseDir/output/VCF_collect"
  input:
  file vcf from vcf2_ch.collect()
  output:
  file "${projectname}_bcf_stats.vchk" into stats_ch
  script:
  """
  bcftools stats -c none ${vcf} > ${projectname}_bcf_stats.vchk
  """
}

// uses HPC modules as cannot get singularity to work this tool
process plots {
  storeDir: "$baseDir/output/VCF_collect"
  input:
  file vchk from stats_ch
  output:
  file "tmp" into move_ch
  script:
  """
  module add bcftools/1.5/gcc
  module add python/anaconda/2.5/2.7
  module add python/anaconda/2019.3/3.7

  plot-vcf -p tmp ${vchk}
  """
}

process move_plots {
  storeDir: "$baseDir/output/VCF_collect/plots"
  input:
  val tmp from move_ch
  output:
  file "venn_bars.indels.png"
  file "venn_bars.snps.png"
  script:
  """
  mkdir $baseDir/output/VCF_collect/plots
  cp $baseDir/output/VCF_collect/tmp/venn_bars.indels.png $baseDir/output/VCF_collect/plots/venn_bars.indels.png
  cp $baseDir/output/VCF_collect/tmp/venn_bars.snps.png $baseDir/output/VCF_collect/plots/venn_bars.snps.png
  rm -fr $baseDir/output/VCF_collect/${tmp}
  rm -fr
  """
}

// add bcftools singularity and targets / regions files for exome
process isec {
  storeDir: "$baseDir/output/VCF_collect/plots"
  input:
  file vcf from vcf3_ch.collect()
  output:
  "${projectname}_merged.vcf.gz" into functotator_ch
  script:
  """
  bcftools isec -c none -O z ${vcf} > ${projectname}_merged.vcf.gz
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
  gatk Funcotator \
   -R /var/spool/mail/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa \
   -V ${vcf} \
   -O ${vcf.baseName}.maf \
   --output-file-format MAF \
   --data-sources-path /var/spool/mail/GATK_functotator_files/funcotator_dataSources.v1.6.20190124g \
   --ref-version hg38
  """
}
