

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
	mkdir -p $baseDir/output/VCF_collect/somalier/somalier_files
	wget -P $baseDir/output/VCF_collect/somalier/somalier_files https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz

	somalier extract -d $baseDir/output/somalier \
	--sites $baseDir/output/VCF_collect/somalier/somalier_files/sites.hg38.vcf.gz \
	-f $genome_fasta \
	${vcf}

	somalier relate --ped $baseDir/chole_batch2.ped  $baseDir/output/somalier/*.somalier

	wget -P $baseDir/output/VCF_collect/somalier/somalier_files \
	https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv

	wget -P $baseDir/output/VCF_collect/somalier/somalier_files \
	https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz?download=1
	tar -xzf $baseDir/output/VCF_collect/somalier/somalier_files/1kg.somalier.tar.gz

	somalier ancestry --labels $baseDir/output/VCF_collect/somalier/somalier_files/ancestry-labels-1kg.tsv \
	$baseDir/output/VCF_collect/somalier/somalier_files/1kg-somalier/*.somalier ++ $baseDir/output/somalier/*.somalier
  """
}
