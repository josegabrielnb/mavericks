#!/usr/local/bin/python3

"""This script uses *.trimmed files of Maverick coordinates to fetch the corresponding
fasta sequences from the NCBI (+25000 nt on each side).

Depends on efetch from the Entrez Programming Utilities"""

import sys, subprocess

my_file = sys.argv[-1]

with open(my_file, 'r') as f:
    lines = f.readlines()
f.closed

data = [line.translate(str.maketrans({'[': '', ']': '', ',': '', '\t': ' ', '\n': ''})) for line in lines]

"""print(*data, sep='\n')"""

for datum in data:
    x = datum.split()
    contig = x.pop(0)
    for i in range(0,int(len(x)/2)):
        coordinate1 = str(int(x[2*i])-25000)
        coordinate2 = str(int(x[2*i+1])+25000)
        submission = 'efetch -db nucleotide -format fasta -id ' + contig + ' ' + '-seq_start ' + coordinate1 + ' ' + '-seq_stop ' + coordinate2
        sequence = subprocess.run(submission, shell = True, stdout = True)
