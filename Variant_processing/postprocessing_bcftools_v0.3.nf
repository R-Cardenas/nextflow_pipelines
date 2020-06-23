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
  storeDir "$baseDir/output/VCF_collect"
  input:
  file vcf from vcf3_ch.collect()
  file index from index3_ch.collect()
  output:
  file "${projectname}_merged.vcf" into split_ch
  script:
  """
  mkdir -p tmp
  bcftools isec -c indels ${vcf} -o ${projectname}_merged.vcf -p tmp
  mv tmp/0001.vcf ./${projectname}_merged.vcf
  rm -fr tmp
  """
}

process split_vcf {
  storeDir "$baseDir/output/VCF_collect/split_vcf"
  input:
  file vcf from split_ch
  output:
  file "*.split.vcf" into (vep_ch, vep_vcf_ch)

  script:
  """
  for file in *.vcf*; do
    for sample in `bcftools query -l \$file`; do
      bcftools view -c1 -s \$sample -o \${file/.vcf*/.\$sample.split.vcf} \$file
    done
  done

  """
}

// you may want to repeat this  but to create the VCF files also
process VEP {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP"
  input:
  file vcf from vep_ch.flatten()
  output:
  file "${vcf.baseName}_VEP.txt" into vep_filter_ch
  script:
  """
  /ensembl-vep/vep -i ${vcf} \
  --dir /var/spool/mail/VEP/.vep \
  -o ${vcf.baseName}_VEP.txt \
  --cache homo_sapiens \
  --sift b \
  --polyphen b \
  --variant_class \
  --af_gnomad \
  --hgvs \
  --domains \
  --tab \
  --show_ref_allele \
  --symbol \
  --nearest gene \
  --verbose
  """
}

// you may want to repeat this  but to create the VCF files also
process VEP2 {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP"
  input:
  file vcf from vep_vcf_ch.flatten()
  output:
  file "${vcf.baseName}_VEP.vcf" into vep2_filter_ch
  script:
  """
  /ensembl-vep/vep -i ${vcf} \
  --dir /var/spool/mail/VEP/.vep \
  -o ${vcf.baseName}_VEP.vcf \
  --cache homo_sapiens \
  --sift b \
  --polyphen b \
  --variant_class \
  --af_gnomad \
  --hgvs \
  --domains \
  --vcf \
  --show_ref_allele \
  --symbol \
  --nearest gene \
  --verbose
  """
}

process vep_filter {
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP"
  input:
  file txt from vep_filter_ch
  output:
  file "${txt.baseName}_filtered.txt" into vep2_ch
  script:
  """
  /ensembl-vep/filter_vep \
  -i ${txt} \
  -o ${txt.baseName}_filtered.txt \
  --format tab \
  --filter "SIFT != tolerated" \
  --filter "SIFT != benign" \
  --filter "SIFT != Tolerated" \
  --filter "SIFT != Benign" \
  --filter "Exome_NFE_AF < 0.1" \
  """
}

process vep_filter2 {
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP"
  input:
  file txt from vep2_filter_ch
  output:
  file "${txt.baseName}_filtered.txt"
  script:
  """
  /ensembl-vep/filter_vep \
  -i ${txt} \
  -o ${txt.baseName}_filtered.txt \
  --format vcf \
  --filter "SIFT != tolerated" \
  --filter "SIFT != benign" \
  --filter "SIFT != Tolerated" \
  --filter "SIFT != Benign" \
  --filter "Exome_NFE_AF < 0.1" \
  """
}



process vep_header {
  storeDir "$baseDir/output/VCF_collect/split_vcf/VEP"
  input:
  file txt from vep2_ch
  output:
  file "${txt.baseName}_noheader.txt"
  script:
  """
  awk '!/\\#/' ${txt} > ${txt.baseName}_noheader.txt
  """
}


// you need to add VCF output
//java -jar /snpEff/snpEff.jar GRCh38.86 chole_batch1_merged.vcf
///gpfs/afm/cg_pipelines/Pipelines/singularity/images/snpeff.simg
