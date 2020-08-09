#!/usr/local/bin/python3

from Bio import Entrez
import time

"""This script downloads Genbank files using the coordinates of intact Mavericks"""

with open('intact_mavericks.txt', 'r') as f:
    mavericks = f.readlines()    
f.closed

Entrez.email = ""

for maverick in mavericks:

    my_seq = maverick.split()[0]

    my_start = maverick.split()[1]

    my_end = maverick.split()[2]

    handle = Entrez.efetch(db="nucleotide", id=my_seq, seq_start=my_start, seq_stop=my_end, rettype="gb", retmode="text")

    with open(my_seq+"_"+my_start+"-"+my_end+".gb", "w") as f:

        print(handle.read(), file=f)
    
    f.close()

    time.sleep(1)


