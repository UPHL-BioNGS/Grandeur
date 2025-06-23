#!/bin/python

##########################################
# written by Erin Young                  #
# for creating summary file for grandeur #
##########################################

import json
import sys

contents=[]

filename = sys.argv[1]
analysis = sys.argv[2]

with open(filename, 'r') as f:
    contents = json.loads(f.read())

sample         = contents['sample']
present        = ','.join(contents['genes']['present'])
absent         = ','.join(contents['genes']['absent'])
susceptibility = ','.join(list(contents['susceptibility'].keys()))

outfile = sample + "_" + analysis + ".tsv"


with open(outfile, 'w') as file:
    file.write("sample\tgenes_present\tgenes_absent\tsusceptibility\n")
    file.write(f"{sample}\t{present}\t{absent}\t{susceptibility}\n")
