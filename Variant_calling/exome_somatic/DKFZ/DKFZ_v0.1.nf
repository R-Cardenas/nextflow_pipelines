/*
 * create a channel for bam files produced by cgpmap_processing pipeline
 */
params.bam = "$baseDir/input/*.cram"
bam_ch = Channel .fromPath( params.bam )

tumor_ch = Channel .fromPath (params.tumor )
normal_ch = Channel .fromPath (params.normal )

println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         GATK Somatic: Mutect2
				 Single samples
				 v0.1
         ===================================



         """
         .stripIndent()

process DKFZ{
  storeDir "$baseDir/ouput/DKFZ"
  input:
  val x from tumor_ch
  val y from normal_ch
  file "${x}.cram" from mutect2_1_ch
  file "${y}.cram" from mutect2_2_ch
  output:

  script:
  """
  run.sh \
    --control-UUID1 --control-bam="${x}.cram" \
    --tumor-UUID2 --tumor-bam="${y}.cram" \
    --delly-svs=/some/path/ICGC_PCA001/delly/bndout/tumor_ICGC_PCA001_vs_control_ICGC_PCA001.bric_embl-delly_2-0-0.20180327.somatic.sv.bedpe.txt \
    --databases=/some/path/dkfz_pipelines_database_files_2018-08-23/ \
    --output=/some/path/ICGC_PCA001/dkfz_pipelines/
  """
}
