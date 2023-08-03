#!/bin/python3

##########################################
# written by Erin Young                  #
# for creating summary files with the    #
# sample id for Grandeur                 #
##########################################

import os
import sys

out = sys.argv[2]
spl = sys.argv[4]

if not os.path.exists(sys.argv[1]):
    print("File " + sys.argv[1] + " does not exist. Exiting.")
    exit()

coms = 0
tabs = 0
with open(sys.argv[1]) as file:
    first_line = file.readline()
    coms = first_line.count('\t')
    tabs = first_line.count('\t')

if tabs > coms:
    delim = '\t'
    print("Predicting tab delimited")
else:
    delim = ','
    print("Predicting comma delimited")

with open(sys.argv[1], 'r') as file:
    lines = file.readlines()
    for line in lines:
        print(line)
        print(line.split(delim))

outfile = open(sys.argv[2], "w")

final_delim = ','
header = 'shouldntexist'

# TODO: turn this into a dict
if sys.argv[3] == 'mlst':
    final_delim = '\t'
    header = 'PubMLST'
    outfile.write('sample\tfilename\tmatching PubMLST scheme\tST\tID1\tID2\tID3\tID4\tID5\tID6\tID7\tID8\tID9\tID10\tID11\tID12\tID13\tID14\tID15\n')
elif sys.argv[3] == 'shigatyper':
    final_delim = '\t'
    header      = 'Number'
elif sys.argv[3] == 'kleborate':
    final_delim = '\t'
    header = 'largest_contig'
elif sys.argv[3] == 'plasmidfinder' :
    final_delim = '\t'
    header = 'Accession number'
elif sys.argv[3] == 'emmtyper':
    outfile.write('sample\tIsolate name\tNumber of BLAST hits\tNumber of clusters\tPredicted emm-type\tPosition(s)\tPossible emm-like alleles\temm-like position(s)\tEMM cluster\n')
    final_delim = '\t'
    header = 'Number of BLAST hits'
elif sys.argv[3] == 'serotypefinder':
    final_delim = '\t'
    header = 'HSP length'

print("Using final delim " + final_delim + " with sample " + spl + " for " + sys.argv[3])

with open(sys.argv[1]) as file:
    lines = file.readlines()
    for line in lines:
        if header in line:
            replace = line.replace(delim, final_delim)
            outfile.write('sample' + final_delim + replace)
        else:
            replace = line.replace(delim, final_delim)
            outfile.write(spl + final_delim + replace)
