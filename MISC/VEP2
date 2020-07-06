
params.vcf= "$baseDir/output/ALL_isec/**/*mutect.vcf"

vcf_ch = Channel. fromPath (params.vcf)


process vcf_stats {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 6
  storeDir "$baseDir/output/ALL_isec/VEP"
  input:
  file vcf from vcf_ch
  output:
  file "${vcf.simpleName}_filtered_AF.vcf" into vcf2_ch
  script:
  """
  bcftools view ${vcf} -i "FORMAT/AF > 0.25 " > ${vcf.simpleName}_filtered_AF.vcf
  """
}


process VEP2 {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/ALL_isec/VEP"
  input:
  file vcf from vcf2_ch
  output:
  file "${vcf.baseName}_VEP.txt" into vep2_ch
  script:
  """
  /ensembl-vep/vep -i ${vcf} \
  --dir /var/spool/mail/VEP_hg19/.vep \
  -o ${vcf.baseName}_VEP.txt \
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


process vep_header {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 6
  storeDir "$baseDir/output/ALL_isec/VEP"
  input:
  file txt from vep2_ch
  output:
  file "${txt.baseName}_noheader.txt"
  script:
  """
  sed -i 's/#Uploaded_variation/Uploaded_variation/g' ${txt}
  awk '!/\\#/' ${txt} > ${txt.baseName}_noheader.txt
  """
}
