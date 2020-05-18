
params.bam ="$baseDir/output/cgpMAP/**/{S5301,S1204,S0108}*.bam"


bam_ch = Channel .fromPath( params.bam )



// the output lookup is wrong
process verifybamid{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/trim/verifyBamID/individual/again"
	input:
	file bam from bam_ch
	output:
	file "*self*"
	file "*depth*"
	file "*log"
	script:
	"""
	samtools index ${bam}
	verifyBamID --vcf $verifybamid --bam ${bam} --out ${bam.simpleName} --maxDepth 1000 --precise --verbose
	"""
}
