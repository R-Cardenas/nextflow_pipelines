[15:39, 10/06/2020] Caz Brown: // Verison v0.1
// Pipeline to process align reads to genome
// source package /tgac/software/testing/bin//nextflow-19.04.1
//Define pathway to fq.gz files
params.reads = "$baseDir/PDPLUG/Outputs/trimgalore/*{1,2}_trim.fq.gz"
//params.reads = "$baseDir/PDPLUG/Outputs/test/*{1,2}_trim.fq.gz"

// Put params into channel
read2_ch = Channel .fromFilePairs (params.reads)

//hisat2 -x /jic/scratch/groups/Christine-Faulkner/Catherine/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Hisat2/hisat2_index -1  -2 test_2_trim.fq.gz -S test.sam

println """\
	\
	\
	\
         ==================================
         RNA-seq - NF
         Mapping and calc counts
         v0.1
         ===================================
         """
         .stripIndent()
//source package…
[15:40, 10/06/2020] Ryan Cardz: wheres the sorted process?
[15:41, 10/06/2020] Caz Brown: omg i deleted it
[15:41, 10/06/2020] Caz Brown: im a fool
[15:41, 10/06/2020] Ryan Cardz: STOP TAPAPING SO FAST!
[15:42, 10/06/2020] Caz Brown: // Verison v0.1
// Pipeline to process align reads to genome
// source package /tgac/software/testing/bin//nextflow-19.04.1
//Define pathway to fq.gz files
params.reads = "$baseDir/PDPLUG/Outputs/trimgalore/*{1,2}_trim.fq.gz"
//params.reads = "$baseDir/PDPLUG/Outputs/test/*{1,2}_trim.fq.gz"

// Put params into channel
read2_ch = Channel .fromFilePairs (params.reads)

//hisat2 -x /jic/scratch/groups/Christine-Faulkner/Catherine/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Hisat2/hisat2_index -1  -2 test_2_trim.fq.gz -S test.sam

println """\
	\
	\
	\
         ==================================
         RNA-seq - NF
         Mapping and calc counts
         v0.1
         ===================================
         """
         .stripIndent()
//source package /nbi/software/testing/bin/hisat-2.0.4
process hisat2{
 storeDir "$baseDir/PDPLUG/Outputs/hisat2"
 input:
 set val(read), file(reads) from read2_ch
 output:
 file "${read}.sam" into fastqc_ch
script:
"""

hisat2  \
-x /jic/scratch/groups/Christine-Faulkner/Catherine/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Hisat2/hisat2_index \
-1 ${reads[0]} \
-2 ${reads[1]} \
-S ${read}.sam
"""
}

process BAM{
  storeDir "$baseDir/PDPLUG/Outputs/BAM"
  input:
  file bam from fastqc_ch
  output:
  file "${bam.simpleName}.bam" into sort_ch
  script:
  """
  source package /nbi/software/production/bin/samtools-1.4.1
  samtools view -S -b ${bam} > ${bam.simpleName}.bam
  """
}

process sort{
  storeDir "$baseDir/PDPLUG/Outputs/BAM"
  input:
  file bam from sort_ch
  output:
  file "${bam.simpleName}_sorted.bam" into counts_ch
  script:
  """
  source package /nbi/software/production/bin/samtools-1.4.1
  samtools sort ${bam} -o ${bam.simpleName}_sorted.bam
  """
}

process stringt{
	storeDir "$baseDir/PDPLUG/Outputs/stringtie"
	input:
	file bam from counts_ch
	output:
	file "${bam.simpleName}.tmp"
	script:
	"""
	source package /nbi/software/testing/bin/stringtie-1.3.3
	stringtie ${bam} -G /jic/scratch/groups/Christine-Faulkner/Catherine/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/TAIR10_GFF3_genes.gff
	"""
}
