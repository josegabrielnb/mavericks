#!/usr/local/bin/python3

"""This script finds blast hits to the Maverick 6 core proteins, must be inside taxon folder
(e.g. Crocodilia/)"""

import sys, os, subprocess

species = sys.argv[-1]

with open(species, 'r') as f:
    spp = f.readlines()
f.closed

for sp in spp:
    folder = sp.rstrip()
    for fasta in os.listdir(folder):
        if fasta.endswith(".fa"):
            print(fasta)
            print("Blasting query proteins against {} ({})...".format(fasta,folder))
            name = fasta.replace(".fa","")
            cmd1 = "tblastn -db " + folder + "/" + fasta + " -query ../Core_genes_query.fa" + " -outfmt '6 qacc pident evalue length' > " + folder + "/" + name + ".blast"
            subprocess.run(cmd1, shell=True)
            print("Done!")
    print()