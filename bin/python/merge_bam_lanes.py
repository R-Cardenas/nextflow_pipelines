# Python script to automate merging of different lanes
## Great toolkit https://github.com/IARCbioinfo/BAM-tricks
# Using info tsv files

import os
import pandas as pd


## import the data sheet provided and extract unique sample names
df = pd.read_csv("../williams_batch2_info.csv")

list1 = list(df['Samples_name'])
unique = list(dict.fromkeys(list1))

print(unique)
# The function that will run the merge bams with the same sample name from CSV
for i in unique:
    script = "samtools merge " + i + "_merged.bam " + i + "*.rename.bam -@ 15"
    print(script)
    #os.system(script)

    # Will rename the SM bam field what is called in the CSV file
    bam_name = i + "_merged.bam"

    script2 = 'samtools view -H ' + bam_name + ' | sed -e "s/SM:[^\\t]*/SM:'+i+'/g" -e "s/ID:[^\\t]*/ID:'+i+'/g" | samtools reheader - ' + bam_name + ' > ' + i + "_merge.rename.bam"
    #os.system(script2)
    print(script2)


