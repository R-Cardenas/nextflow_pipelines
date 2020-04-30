/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/output2/cgpMAP/*.rename.bam"
params.index = "$baseDir/output2/cgpMAP/*.bam.bai"
params.bas = "$baseDir/output2/cgpMAP/*.bam.bas"
bam_ch = Channel .fromPath( params.bam )
index_ch = Channel .fromPath( params.index )
bas_ch = Channel .fromPath( params.bas )

bam_ch.into { bam2_ch; bam3_ch }
index_ch.into { index2_ch; index3_ch }
bas_ch.into { bas2_ch; bas3_ch }

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
				 v0.1
         ===================================



         """
         .stripIndent()


process cgpwxs2 {
	storeDir "$baseDir/output/cgpwxs"
	stageInMode 'copy'

	input:
	val tumor from tumor_ch
	val normal from normal_ch
	file "${tumor}.bam" from bam2_ch
	file "${tumor}.bam.bai" from index2_ch
	file "${tumor}.bam.bas" from bas2_ch
	file "${normal}.bam" from bam3_ch
	file "${normal}.bam.bai" from index3_ch
	file "${normal}.bam.bas" from bas3_ch
	output:
	file "WXS_${tumor}_vs_${normal}.result.tar.gz"
	script:
	"""
	singularity exec --cleanenv \
	--home $baseDir/output/cgpwxs \
	--bind /gpfs/afm/cg_pipelines/Pipelines/singularity/genomes:/var/spool/ref:ro \
	--bind $baseDir/output2/cgpMAP:/var/spool/data:ro \
	/gpfs/afm/cg_pipelines/Pipelines/singularity/images/cgpwxs_3.1.6.img \
	ds-cgpwxs.pl \
	 -reference /var/spool/ref/cgpwgs_ref/GRCh38/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
   -annot /var/spool/ref/cgpwgs_ref/GRCh38/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz \
	 -snv_indel /var/spool/ref/cgpwgs_ref/GRCh38/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz \
	 -tumour ${tumor}.bam \
	 -tidx ${tumor}.bam.bai \
	 -normal ${normal}.bam \
	 -nidx ${normal}.bam.bai \
	 -exclude NC_007605,hs37d5,GL% \
	"""
}
