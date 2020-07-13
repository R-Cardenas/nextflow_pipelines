
params.vcf= "$baseDir/output/ALL_isec_indels/*vcf"

vcf_ch = Channel. fromPath (params.vcf)

process awks_stats {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/ALL_isec_indels"
  input:
  file vcf from vcf_ch
  output:
  file "${vcf.simpleName}_indels.vcf" into vep_ch
  script:
  """
  bcftools view -v indels ${vcf} > ${vcf.simpleName}_indels.vcf
  """
}

process VEP{
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/ALL_isec_indels"
  input:
  file vcf from vep_ch
  output:
  file "${vcf.simpleName}_VEP.txt"
  script:
  """
  /ensembl-vep/vep  -i ${vcf} \
    -o ${vcf.simpleName}_VEP.txt \
    --cache homo_sapiens \
    --dir_cache /var/spool/mail/VEP_hg19/.vep \
    --offline \
    --fasta $VEP_fasta \
    --fork 5 \
    --cache homo_sapiens \
    --sift p \
    --polyphen p \
    --variant_class \
    --af_gnomad \
    --no_stats \
    --tab \
    --show_ref_allele \
    --symbol \
    --verbose \
    --domains \
    --regulatory
  """
}
