

/*
 * create a channel for fastq pairs
 */
params.bam = "/gpfs/home/aft19qdu/scratch/test/output/cgpMAP/**/*.bam"

Channel
	.fromPath( params.bam )
	.collect ()
	.set {bam_ch, bam2_ch}


println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         Freebayes - Germline
         v0.2
         ===================================



         """
         .stripIndent()


process somalier{
  storeDir "$baseDir/output/somalier"
  input:
  file bam from bam_ch.collect()
  output:
	file "*.html"
  script:
  """
	for f in *.bam; do
    somalier extract -d $baseDir/output/somalier/extracted/ --sites $baseDir/sites.vcf.gz -f $genome_fasta \$f
	done

	somalier relate --ped $baseDir/chole_batch2.ped  $baseDir/output/somalier/extracted/*.somalier
  """
}

process ancestry{
  storeDir "$baseDir/output/somalier"
  input:
  file bam from bam2_ch.collect()
  output:
	file "*.html"
  script:
  """
	somalier ancestry --labels $baseDir/ancestry-labels-1kg.tsv \
	$baseDir/1kg-somalier/*.somalier ++ $baseDir/output/somalier/extracted/*.somalier
  """
}
