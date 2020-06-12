
params.vcf = "$baseDir/output/VCF_collect2/VCF_collect/plots/fam_vcf_maf/intersected_by_fam/*.vcf.gz"

vcf_ch = Channel. fromPath (params.vcf)


process VEP{
  storeDir "$baseDir/output/VCF_collect2/VCF_collect/plots/fam_vcf_maf/VEP"
  input:
  file vcf from vcf_ch
  output:
  file "${vcf.simpleName}_VEP.txt"
  script:
  """
  vep -i ${vcf} \
    -o ${vcf.simpleName}_VEP.txt \
    --cache homo_sapiens \
    --dir_cache /var/spool/mail/VEP/.vep \
    --polyphen p \
    --sift p \
    --tab \
    --show_ref_allele \
    --custom /var/spool/mail/VEP/gnomad_hg38/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_NFE \
    --hgvs \
    --domains \
    --symbol \
    --coding_only \
    --verbose
  """
}
