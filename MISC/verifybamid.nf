
params.bam ="$baseDir/output/BAM/merged_lanes/*.rmd.bam"
params.idx ="$baseDir/output/BAM/merged_lanes/*.bai"

bam_ch = Channel .fromPath( params.bam )
idx_ch = Channel .fromPath( params.idx )



process verifybamid{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/trim/verifyBamID"
	input:
	file bam from bam_ch
	file idx from idx_ch.collect()
	output:
	file "verifybam*"
	script:
	"""
	verifyBamID --vcf $verifybamid --bam ${bam} --out ${bam.simpleName} --maxDepth 1000 --precise --verbose
	"""
}
