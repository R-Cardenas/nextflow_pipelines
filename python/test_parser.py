#!/usr/bin/env python3
# Python script to run gatk CombineGVCFs
import argparse
import os

#Create the input called scquare
parser = argparse.ArgumentParser()
parser.add_argument("-V", dest="square", type=str,
                    help="display a square of a given number")
parser.add_argument("-O", dest="triangle", type=str,
                    help="display a triangle to output")

args = parser.parse_args()

#Append square agruement to example obj
example = args.square
output = args.triangle

x = example.split()

# For loop adding samples with a V and new line
y = ""

for i in x:
    y += "-V " + i + " "

#Put this together with the GATK tool to run run the script we want
script = "gatk CombineGVCFs -R /var/spool/mail/hg38/UCSC/WholeGenomeFasta/genome.fa " + y + "-O " + output
script = str(script)

print(script)
#Run this scrip by the OS
os.system(script)
