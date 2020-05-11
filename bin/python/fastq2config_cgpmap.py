#!/usr/bin/env python

import argparse

# # set args parse to allow input from terminal

parser = argparse.ArgumentParser()

parser.add_argument('--fq1', required=True)
parser.add_argument('--fq2', required=True)
parser.add_argument('--n1', required=True)
parser.add_argument('--n2', required=True)
parser.add_argument('--o', required=True)

args = parser.parse_args()
f1 = args.fq1
f2 = args.fq2
n1 = args.n1
n2 = args.n2
o1 = args.o

# open the files
file1 = open(str(f1), "r")
file1 =str(file1.read())
print(file1)

file2 = open(str(f2), "r")
file2 =str(file2.read())
print(file2)

# remove white space and split string by ':'
file1 = str.strip(file1)
file2 = str.strip(file2)
f1_list = str.split(file1,":")
f2_list = str.split(file2,":")


## Setting paramters for the config file - see for more details https://github.com/cancerit/PCAP-core/wiki/File-Formats-groupinfo.yaml

# Sample names from both files and ensure they are same for both
name1 = str.split(f1,"_")
name1 = name1[0]
name2 = str.split(f2,"_")
name2 = name2[0]

if name1 == name2:
    print("names are identical")
else:
    raise SyntaxError('Fastq files have different names')
    exit
# ## Set other parameters from input within txt file
SM = "SM: " + name1 + "\n"
PL = "    PL: ILLUMINA" + "\n"
RG = "READGRPS:" + "\n"
#FQ1
FQ1 = "  " + n1 + ":\n"
LB1 = "    LB: " + name1 + "_" + f1_list[9] + "\n"
PU1 = "    PU: " + f1_list[2] + "." + f1_list[3] + "\n"
#FQ2
FQ2 = "  " + n2 + ":\n"
LB2 = "    LB: " + name2 + "_" + f2_list[9] + "\n"
PU2 = "    PU: " + f2_list[2] + "." + f2_list[3] + "\n"
# ## Forming the final YAML file for cgpmap and save file
yaml = SM + RG + FQ1 + PL + LB1 + PU1 + FQ2 + PL + LB2 + PU2
# print(yaml)
filename  = o1
f= open(filename,"w+")
f.write(yaml)
# ## Rationale for using seq primer and sample id for LB input
### MarkDuplicates uses the LB field to determine which read groups might contain molecular duplicates,
# in case the same DNA library was sequenced on multiple lanes.
# Therefore if the same library will have same sequencing primers..EZ
# ## regex for pulling index
# grep ":[C,A,T,G]*[+][C,A,T,G]" new.txt | head -1
# fqtools -d header S0102-P_EKDN200000457-1A_HYFMTDSXX_L1_1.fq.gz | grep ":[C,A,T,G]*[+][C,A,T,G]" | head -1
