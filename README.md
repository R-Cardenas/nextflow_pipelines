# Variant call pipeline (nextflow)

## Overview

This github houses the nextflow scripts required to process DNA-seq (WGS and exome) samples from fastq files to calling variants. The pipelines are formed of different modules that can be used in different combinations to fulfil the requried analysis. Below shows a picture-graphic of the modules available currently.

![figure-1](images/exome_flow_summary.png)


As depicted in figure-1, each module is connected through the last file is produced. For example, the final files that are produced by cgpMAP are bam files with the file extension ''*.merged.bam'. When run, the connecting modules search for files ending in '*.merged.bam' and then start. (Currently) Once the pipeline has determined the number of input files and started, none more can be added. This means that the modules are required to run in succession, which reduces the speed of the overall pipeline.

Each module is comprised of various QC steps to ensure high quality processing of samples. The following sections describes the steps within each module and how it is expected to work.


## Mapping-exome module

The following module performs QC on fastq files and maps them using the sanger cgpMAP container (bwa-mem). Subsequent BAM files are de-duplicated and QC'ed (insert size, hybrid stats, alignment stats). BAM files are collected in order to merge the lanes.

### Trim-galore

Trimgalore requires the index primers to be supplied in order to trim fastQ files. These sequences can be supplied in the input.config file. FastQC option is enabled in order to produce report.

### CgpMAP and fqtools





For the following section there are several files that are required to be supplied:

  - exome target/bait intervals (hybrid stats)

![figure-2](images/mapping_exome.png)
