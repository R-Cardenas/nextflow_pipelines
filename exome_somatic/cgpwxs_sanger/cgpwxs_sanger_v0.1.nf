/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output/hg38_decoy/aligned_sorted/*.rename.bam"
params.index = "$baseDir/output/hg38_decoy/aligned_sorted/*.rename.bam.bai"
bam_ch = Channel .fromPath( params.bam ) .ifEmpty('No bam files detected') .println()
index_ch = Channel .fromPath( params.index ) .ifEmpty('No index files detected') .println()

// Split channels so can be used twice
bam_ch.into { bam2_ch; bam3_ch }
index_ch.into { index2_ch; index3_ch }

// Define Tumor vs Normal Variables
tumor_ch = Channel .from (params.tumor ) .ifEmpty('Please specify tumor samples') .println()
normal_ch = Channel .from (params.normal ) .ifEmpty('Please specify normals samples') .println()

println """\
	\
	\
	\
         ==================================
         E X O M E - N F   S O M A T I C
         cgpwxs_3.1.6.img
				 Tumour matched
				 v0.1
         ===================================



         """
         .stripIndent()

process cgpwxs2 {
	storeDir "$baseDir/output/hg38_decoy/cgpwxs"
	input:
	val tumor from tumor_ch
	val normal from normal_ch
	file "${tumor}.rename.bam" from bam2_ch
	file "${tumor}.rename.bam.bai" from index2_ch
	file "${normal}.rename.bam" from bam3_ch
	file "${normal}.rename.bam.bai" from index3_ch
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
	 -tidx $baseDir/output/hg38_decoy/aligned_sorted/${tumor}.rename.bam.bai \
	 -normal $baseDir/output/hg38_decoy/aligned_sorted/${normal}.rename.bam \
	 -nidx $baseDir/output/hg38_decoy/aligned_sorted/${normal}.rename.bam.bai \
	 -exclude NC_007605,hs37d5,GL% \
	 -outdir $baseDir/output/hg38_decoy/cgpwxs
	"""
}
