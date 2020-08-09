#!/usr/local/bin/python3

"""This script builds HMMs of the first 200-bp of the 5' inverted repeat of Maverick alignments"""

import sys, os, subprocess, re

with open("coelacanthi_species.txt", 'r') as f:
            spp = f.readlines()
f.closed

for sp in spp:

    species = sp.strip()

    print("Building HMM for {}...".format(species))

    os.chdir(species)

    cmd1 = 'cat *.fa > ' + species + '_mavericks.fa'
    cmd2 = 'mafft --thread 4 ' + species + '_mavericks.fa > ' + species + '_mavericks.mafft.aln'
    cmd3 = 'extractalign -sequence ' + species + '_mavericks.mafft.aln -regions 1-200 -outseq ' + species + '_5prime200.aln'
    cmd4 = 'hmmbuild ' + species + '_5prime200.hmm ' + species + '_5prime200.aln'

    print("Concatenating sequences")
    subprocess.run(cmd1,shell=True)
    print("Aligning sequences")
    subprocess.run(cmd2,shell=True)
    print("Extracting sequences")
    subprocess.run(cmd3,shell=True)
    print("Building HMM")
    subprocess.run(cmd4,shell=True)
    print("Done!\n")

    os.chdir('..')

