# Python script to automate merging of different lanes
## Great toolkit https://github.com/IARCbioinfo/BAM-tricks
# Using info tsv files

import os
import pandas as pd
from pathlib import Path
import time

## import the data sheet provided and extract unique sample names
df = pd.read_csv("williams_batch2_info.csv")

list1 = list(df['Samples_name'])
unique = list(dict.fromkeys(list1))

print(unique)
# The function that will run the merge bams with the same sample name from CSV
for i in unique:

    script = '"samtools merge ' + i + '_merged.bam ' + i + '*.sorted.bam -@ 5"'
    script2 = 'bsub -q long-ib -M 20000 -R"rusage[mem=20000]" -J MERGE ' + script
    print(script2)
    os.system(script2)
    time.sleep(180)

def loop():
    for i in unique:
        file = i + "_merged.bam"
        if os.path.isfile(file):
            print ("File exist")
        else:
            print('file does not exist yet')
            time.sleep(360)
            loop()

def size():
    for i in unique:
        file = i + "_merged.bam"
        size1 = os.stat(file).st_size
        time.sleep(180)
        size2 = os.stat(file).st_size
        if size1 == size2:
            print("files are the same size")
            os.system("echo 'size1, size2 " + size1 + size2 + file + "' >> size.log")
        else:
            size()

loop()
size()
time.sleep(500)
