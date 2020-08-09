#!/bin/local/python3

"""Get cDNAs in fasta format from GeneWise files"""

import os, re

with open('Genewise_ATP_cDNAs.output', 'r') as f:
    lines = f.readlines()
f.close()

index = 0
name = ''

for line in lines:

    if re.match("Target Sequence", line):

        name = line[21:]
        
        name = ">" + name.replace(".fasta", "")
        
        print(name.strip())

    if re.match(r'.*Alignment.*', line):

        index = 1

    elif re.match(r'\/\/', line):

        index = 0

        name = ""
    
    if index == 1:

        seqs = line[20:]

        seqs = seqs.strip().replace(" ", "")

        if seqs.islower() or re.match(r'.*Intron.*', seqs):

            seqs = re.sub('<', '', seqs)
            seqs = re.sub('>', '', seqs)
            seqs = re.sub('-', '', seqs)
            seqs = re.sub('\[', '', seqs)
            seqs = re.sub('\]', '', seqs)
            seqs = re.sub('Intron', '', seqs)
            seqs = re.sub(r'[0-9]*', '', seqs)
            seqs = re.sub(r'[ACGT]*', '', seqs)

            print(seqs)