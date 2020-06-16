

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
	mkdir -p bin
	wget -P bin/ https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz

	for f in *.cram; do
    somalier extract -d extracted/ --sites bin/sites.hg38.vcf.gz -f -f $genome_fasta \$f
	done

	somalier relate --ped $baseDir/bin/chole_batch2.ped  bin/*.somalier

	wget -P ancestry_files https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv

	wget https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz
	tar -xzf 1kg.somalier.tar.gz

	somalier ancestry --labels ancestry_files/ancestry-labels-1kg.tsv \
	1kg-somalier/*.somalier ++ extracted/*.somalier
  """
}
