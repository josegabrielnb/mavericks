#!/usr/local/bin/python3

"""Perform remote tblastn searches against the nucleotide non-redundant database"""

import subprocess, time

for i in range(1,43):

    time.sleep(5)

    cmd = 'tblastn -db nr -query ' + str(i) + '.fa -evalue 1e-100 -outfmt "6 qacc sacc evalue sseq" -entrez_query "vertebrates[Organism]" -remote -out pPOLB_tblastn_' + str(i) + '.txt'

    subprocess.run(cmd, shell = True)
    
