#!/usr/bin/env python

import argparse

# positional args

parser = argparse.ArgumentParser()

parser.add_argument('--fq1', required=True)
parser.add_argument('--fq2', required=True)

args = parser.parse_args()
f1 = args.fq1
f2 = args.fq2
