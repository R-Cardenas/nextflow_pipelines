/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */

params.bam = "$baseDir/output/hg38_decoy/aligned_sorted/*.rename.{bam.bai}"

bam_ch = Channel .fromFilePairs( params.bam )


// Split channels so can be used twice
bam_ch.into { bam2_ch; bam3_ch }

// Define Tumor vs Normal Variables
tumor_ch = Channel .from (params.tumor )
normal_ch = Channel .from (params.normal )

println """\
	\
	\
	\
         ==================================
         E X O M E - N F   S O M A T I C
         cgpwxs_3.1.6.img
				 Tumour matched
				 v0.2
         ===================================



         """
         .stripIndent()

process cgpwxs2 {
	storeDir "$baseDir/output/hg38_decoy/cgpwxs"
	input:
	val tumor from tumor_ch
	val normal from normal_ch
 	set val(xx), file("${tumor}.rename.bam"), file("${tumor}.rename.bai") from bam2
	set val(yy), file("${normal}.rename.bam"), file("${normal}.rename.bai") from bam3
	output:
	file "WXS_${tumor}_vs_${normal}.result.tar.gz"
	script:
	"""
	singularity exec --cleanenv \
	--home $baseDir \
	--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
	--bind $baseDir/output/cgpMAP:/var/spool/data:ro \
	/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwxs_3.1.6.img \
	ds-cgpwxs.pl \
	 -reference /var/spool/ref/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
   -annot /var/spool/ref/cgpwgs_ref/GRCh38/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz \
	 -snv_indel /var/spool/ref/cgpwgs_ref/GRCh38/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz \
	 -tumour $baseDir/output/hg38_decoy/aligned_sorted/${tumor}.rename.bam \
	 -tidx $baseDir/output/hg38_decoy/aligned_sorted/${tumor}.rename.bai \
	 -normal $baseDir/output/hg38_decoy/aligned_sorted/${normal}.rename.bam \
	 -nidx $baseDir/output/hg38_decoy/aligned_sorted/${normal}.rename.bai \
	 -exclude NC_007605,hs37d5,GL% \
	 -outdir $baseDir/output/hg38_decoy/cgpwxs
	"""
}
