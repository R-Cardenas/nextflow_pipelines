# Python script to automate merging of different lanes
## Great toolkit https://github.com/IARCbioinfo/BAM-tricks
# Using info tsv files

import os
import pandas as pd


## import the data sheet provided and extract unique sample names
df = pd.read_csv("williams_batch2_info.csv")

list1 = list(df['Samples_name'])
unique = list(dict.fromkeys(list1))

print(unique)
# The function that will run the merge bams with the same sample name from CSV
for i in unique:
    script = "samtools merge " + i + "_merged.bam " + i + "*.rmd.bam -@ 15"
    print(script)
    os.system(script)
