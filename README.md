# Variant call pipeline (nextflow)

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [Overview](#Overview)
* [Mapping-exome module](#Mapping-exome)
* [Freebayes individual](#Freebayes)
* [GATK Cohort](#GATK)


## Overview

This github houses the nextflow scripts required to process DNA-seq (WGS and exome) samples from fastq files to calling variants. The pipelines are formed of different modules that can be used in different combinations to fulfil the requried analysis. Below shows a picture-graphic of the modules available currently.

![figure-1](images/exome_flow_summary.png)


As depicted in figure-1, each module is connected through the last file is produced. For example, the final files that are produced by cgpMAP are bam files with the file extension ''*.merged.bam'. When run, the connecting modules search for files ending in '*.merged.bam' and then start. (Currently) Once the pipeline has determined the number of input files and started, none more can be added. This means that the modules are required to run in succession, which reduces the speed of the overall pipeline.

Each module is comprised of various QC steps to ensure high quality processing of samples. The following sections describes the steps within each module and how it is expected to work.


## Mapping-exome

Link to current cgpMAP nextflow script [here](cgpmap/cgpmap_preprocessing_v0.7.nf)

The following module performs QC on fastq files and maps them using the sanger cgpMAP container (bwa-mem). Subsequent BAM files are de-duplicated and QC'ed (insert size, hybrid stats, alignment stats). BAM files are collected in order to merge the lanes. Following is a brief breakdown of the processes within the pipeline.

![figure-2](images/mapping_exome.png)

### Trim-galore

Trimgalore requires the index primers to be supplied in order to trim fastQ files. These sequences can be supplied in the input.config file. FastQC option is enabled in order to produce report.

### CgpMAP and fqtools

CgpMAP does not annotate the bam files with information from the fastQ header. Therefore following trimming the fastq files are passed to FQtools which will extract the header for each read pair and a python script (fastq2config_cgpmap.py) to write this into a YAML file. The YAML is passed to cgpMAP which allows the correct headers to be assigned. An example of a YAML output is shown below:

```
SM: S1001
READGRPS:
  S1001_EKDN200000467-1A_HYFMTDSXX_L1_1.fq.gz:
    PL: ILLUMINA
    LB: S1001_TTATCGGC+GATCATCC
    PU: HYFMTDSXX.1
  S1001_EKDN200000467-1A_HYFMTDSXX_L1_2.fq.gz:
    PL: ILLUMINA
    LB: S1001_TTATCGGC+GATCATCC
    PU: HYFMTDSXX.1
```

### BAM QC

The following tools were used:
  - Picard remove duplicates
  - Picard insert size
  - picard hybrid stats
  - picard bam alignment stats

### Merge lanes

Following the removing of duplicates from BAM files, the lanes are merged with supplied information from an excel sheet that is processed by [python script](bin/python/merge_bam_lanes_2.py). An example of the excel sheet layout is shown [here](MAP/cgpmap/williams_batch2_info.csv). For future work it would be easier to merge all bam files that have the same sample name.

### Exome mapping problems and workarounds

While constructing this pipeline there were a number of issues that were resovled using hacky techniques. The first is not confined to this module, but in fact is a pipeline wide issue. Error 255/130:
```
{put in error}
```
This was a particular issue with singularity images harboring picard and samtools (theyre together) and cgpMAP. For samtools/picard the issue appeared to be an issue of resource - ie too many processes accessing the image at once resulted in failure. This was a similar issue in cgpMAP, however the problem was cgpMAP accessing the reference genomes to unzip. If too many processes did this at once resulted in error.

In order to circumvent this, job submissions were limited to 10 at any one time (per process)( [nextflow.config:](config/nextflow_v0.3.config): queueSize = 10) and to space submissions out ([nextflow.config:](config/nextflow_v0.3.config) submitRateLimit = '1 / 1min'). This reduced the number of errors, though an error error strategy is still required for these processes to retry:

```
errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
maxRetries 6
```

This error strategy allows task attempts to have varied times when they are retried to increase the chances of processes having resource available when resubmitted.


### To be done

- script to download the scripts required for all pipeline before pipeline starts (easy github download)
- clean up channel names.. messy because of extensive changing
- merge lanes to be optional.. not all pipelines require merging of lanes.. a seperate nextflow script with it excluded will be the easiest.
- verifyBamID(1) to be implemented. The tool is now working just needs to be added.


## Freebayes_(individual)

CgpMap produces *.merge.bam files within $baseDir/output/merge_lanes directory that are detected by the Freebayes nextflow pipeline.

Each bam file is run independantly by Freebayes and filtered using the freebayes recommended parameters:


```
vcffilter -f \
  "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
  ${vcf} > ${vcf}_filtered_freebayes.vcf
```

VCF files are merged into one and moved to the VCF_collect folder, awaiting variant post processing (!!LINK)

Problems/ to be done: NONE

## GATK(cohort_mode)

CgpMap produces *.merge.bam files within $baseDir/output/merge_lanes directory that are detected by the GATK nextflow pipeline.
