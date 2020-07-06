/*
 * create a channel for fastq pairss
 */

// Input Reads
params.read1 = "/gpfs/afm/cg_pipelines/Pipelines/Chole_batch1_FASTQ/*{1,2}.fq.gz"


read1_ch = Channel .fromFilePairs( params.read1 )
read1_ch.into { read2_ch; read3_ch }

params.csv = "$baseDir/bin/williams_batch2_info.csv"
csv_ch = Channel .fromPath( params.csv )

println """\
	\
	\
	\
         ==================================
         E X O M E - N F   P I P E L I N E
         Mapping and Bam Processing
         v0.3
         ===================================



         """
         .stripIndent()


myLongCmdline = "git clone https://github.com/R-Cardenas/nextflow_pipelines.git"
result = myLongCmdline.execute().text

// link for adding fastp if the pipe works https://github.com/nextflow-io/nextflow/issues/682
process trim_galore{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	stageInMode = 'copy' // trim_galore doesnt like sym/hardlinks..
  storeDir "$baseDir/output/cgpMAP/trim_galore"
	input:
	tuple val(read2), file(reads) from read2_ch
	output:
	file "${reads[0].simpleName}.fq.gz" into (read5_ch, read7_ch)
	file "${reads[1].simpleName}.fq.gz" into (read10_ch, read12_ch)
	script:
	"""
	mkdir -p $baseDir/logs

	trim_galore --paired \
	--fastqc --illumina \
	${reads[0]} ${reads[1]}

	mv ${reads[0].simpleName}_val_1.fq.gz ${reads[0].simpleName}.fq.gz
	mv ${reads[1].simpleName}_val_2.fq.gz ${reads[1].simpleName}.fq.gz

	rm -fr ${reads[0]} # remove the copied files to prevent memory loss
	rm -fr ${reads[1]}
	"""
}


process fqtools{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/cgpMAP/trim_galore"
	input:
	file read1 from read7_ch
	file read2 from read12_ch
	output:
	file "${read1}.yaml" into yaml_ch
	file("fqtools_WARNING_?.txt") optional true
	script:
	"""
	fqtools -d header ${read1} | grep ":[C,A,T,G]*[+][C,A,T,G]" | head -1 > ${read1.simpleName}.txt

	### Counts lines in file and will repeat if it is empty
	words=`wc -l ${read1.simpleName}.txt  | awk '{print \$1}'`;
	if [ \$words -eq 0 ]
	then
	fqtools -d header ${read1} | head -1 > ${read1.simpleName}.txt
	else
	echo 'alls good!'
	fi

	### Counts lines in file and will repeat if it is empty
	fqtools -d header ${read2} | grep ":[`§`,A,T,G]*[+][C,A,T,G]" | head -1 > ${read2.simpleName}.txt
	words=`wc -l ${read2.simpleName}.txt  | awk '{print \$1}'`;
	if [ \$words -eq 0 ]
	then
	fqtools -d header ${read2} | head -1 > ${read2.simpleName}.txt
	else
	echo 'alls good!'
	fi

	cp fqtools_WARNING_?.txt $baseDir/logs
	python $baseDir/nextflow_pipelines/bin/python/fastq2config_cgpmap.py \
	--fq1 ${read1.simpleName}.txt --fq2 ${read2.simpleName}.txt \
	--n1 ${read1} --n2 ${read2} --o ${read1}.yaml
	"""
}

process cgpMAP {
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/cgpMAP/${read1.simpleName}"
  input:
	val read1 from read5_ch
	val read2 from read10_ch
	val yaml from yaml_ch.collect()
  output:
  file "*.bam" into cgp_ch
  script:
  """

  name=\$(echo '${read2}' | sed -e 's/.*[/]//' -e 's/_.*//')

  ds-cgpmap.pl  \
  -outdir $baseDir/output/cgpMAP/${read1.simpleName} \
  -r $cgpmap_genome \
  -i $cgpmap_index \
  -s \$name \
  -t 5 \
	-g ${read1}.yaml \
  ${read1} ${read2}

	mv $baseDir/output/cgpMAP/${read1.simpleName}/*.bam \
	$baseDir/output/cgpMAP/${read1.simpleName}/${read1.simpleName}.bam

	echo 'fq1: ${read1} fq2: ${read2} bam_name: ${read1.simpleName}' >> $baseDir/${projectname}_cgpmap_samples.log
	ls -l  ${read1} >> $baseDir/logs/symbolic_test_fastq.log
	ls -l  ${read2} >> $baseDir/logs/symbolic_test_fastq.log
  """
}


process sam_sort {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
  storeDir "$baseDir/output/BAM/sorted"
  input:
  file bam from cgp_ch
  output:
  file "${bam.simpleName}.sorted.bam" into dup_ch
  script:
  """
	mkdir -p tmp
  picard SortSam I=${bam} O=${bam.simpleName}.sorted.bam SORT_ORDER=coordinate TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
  """
}


process picard_pcr_removal {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/BAM/merged_lanes"
  input:
  file bam from dup_ch.flatten()
  output:
  file "${bam.simpleName}.rmd.bam" into (index1_ch, index_2ch, hs_ch, bam10_ch, bam11_ch, bam12_ch)
	file "${bam.simpleName}.log"
  script:
  """
	mkdir -p tmp
  picard MarkDuplicates I=${bam} O=${bam.simpleName}.rmd.bam M=${bam.simpleName}.log TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process bam_index {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/BAM/merged_lanes"
  input:
  file bam from index1_ch
  output:
  file "${bam}.bai" into index_3ch

  script:
  """
	mkdir -p tmp
  picard BuildBamIndex \
	I=${bam} \
	O=${bam}.bai \
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process collect_insert_size {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/BAM/insert_size"
  input:
  file bam from index_2ch
  output:
  file "${bam.simpleName}_insert_size.txt"
	file "${bam.simpleName}_insert_size.txt"
  script:
  """
	mkdir -p tmp
  picard CollectInsertSizeMetrics I=${bam} H=${bam.simpleName}_histogram.pdf \
  O=${bam.simpleName}_insert_size.txt M=0.5 TMP_DIR=tmp
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process hybrid_stats {
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/BAM/hybrid_stats"
  input:
  file bam from hs_ch
  output:
  file "${bam.simpleName}_hs_metrics.txt"
  script:
  """
	mkdir -p tmp
  picard CollectHsMetrics I=${bam} O=${bam.simpleName}_hs_metrics.txt \
  R=$genome_fasta \
  BAIT_INTERVALS=$bait_interval \
  TARGET_INTERVALS=$target_interval \
	TMP_DIR=tmp
	rm -fr tmp
  """
}

process alignment_stats{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/BAM/alignment_stats"
	input:
	file bam from bam10_ch
	output:
	file "${bam.simpleName}_align_stats.txt" into verify_ch
	script:
	"""
	mkdir -p tmp
	picard CollectAlignmentSummaryMetrics \
  R=$genome_fasta \
	I=${bam} \
	O=${bam.simpleName}_align_stats.txt \
	TMP_DIR=tmp
	rm -fr tmp
	"""
}

process verifybamid{
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 6
	storeDir "$baseDir/output/BAM/verifyBamID"
	input:
	file bam from bam11_ch
	file idx from index_3ch.collect()
	output:
	file "*.depthRG"
	file "*.depthSM"
	file "*.log"
	file "*.selfRG"
	file "*.selfSM"
	script:
	"""
	verifyBamID --vcf $verifybamid --bam ${bam} --out ${bam.simpleName} --maxDepth 1000 --precise --verbose
	"""
}



process somalier{
  storeDir "$baseDir/output/BAM/somalier"
  input:
  file bam from bam12_ch.collect()
  output:
	file "*.html"
  script:
  """
	mkdir -p bin
	wget -P bin/ https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz

	for f in *.bam; do
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



workflow.onComplete {
	// create log files and record output
	process finish{
		script:
		""" echo ' Pipeline cgpmap v0.3 completed
		Project: $projectname
		Time: ${nextflow.timestamp}

	 trimmomatic - completed
	 (adapter1: $adapter1
		adapter2: $adapter2)
	 cgpmap - completed
	 (reference = $cgpmap_genome)
	 (index = $cgpmap_index)
	 samsort - completed
	 picard pcr removal - completed
	 picard rename bam - completed
	 samtools bam index - completed
	 picard collect insert size - completed
	 picard collect HS stats - completed
	 (target_intervals = $target_interval)
	 (bait_intervals = $target_interval)
	 merge bams - completed
	 ' >> $baseDir/${projectname}_log.txt

	 mail -s "cgpMAP successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_log.txt
	 	mail -s "cgpMAP info" aft19qdu@uea.ac.uk < $baseDir/${projectname}_cgpmap_samples.log
	 """
	}
}

workflow.onError {
	process finish_error{
		script:
		"""
		echo 'Pipeline cgpmap v0.3 FAILED
		Project: $projectname
		Time: ${nextflow.timestamp}

		Error:
		${workflow.errorMessage}' >> $baseDir/${projectname}_error.txt

	  mail -s "cgpMAP successful" aft19qdu@uea.ac.uk < $baseDir/${projectname}_error.txt
	  """
	}
}
