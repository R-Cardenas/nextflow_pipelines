# Python script to automate merging of different lanes
## Great toolkit https://github.com/IARCbioinfo/BAM-tricks
# Using info tsv files

import os
import pandas as pd
import glob


## import the data sheet provided and extract unique sample names
df = pd.read_csv("../williams_batch2_info.csv")

for ind in df.index:
    file = df['FASTQ_names'][ind] + '*'
    file2 = glob.glob(file)
    file2 = str(file2[0])
    script2= 'mkdir -p tmp'
    #os.system(script2)
    script3 = "picard AddOrReplaceReadGroups I=" + file2 + \
        " RGID=" + str(df['Read_ID'][ind]) + '.' + str(df['Lane'][ind]) + \
        " RGLB=lib1" + " RGPL=illumina" + \
        " RGSM=" + str(df['Samples_name'][ind]) + \
        " O=" + str(df['FASTQ_names'][ind]) + ".rename.bam " + \
        "TMP_DIR=tmp" + " RGPU=twat"
    print(script3)
    #os.system(script3)




