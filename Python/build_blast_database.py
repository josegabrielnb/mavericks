#!/usr/local/bin/python3

"""This script aims to find Maverick TIRs in candidate regions of the genome.
Blast command line executables must be installed system wide."""

import sys, os, subprocess

folder = sys.argv[-1]

files = os.listdir(folder)
#files.remove(".DS_Store")

for fasta in files:
    cmd1 = "makeblastdb -in " + folder + "/" + fasta + " -dbtype nucl -parse_seqids"
    subprocess.run(cmd1, shell=True)

for fasta in files:
    name = fasta.replace(".fa","")
    cmd2 = "blastn -db " + folder + "/" + fasta + " -query " + folder + "/" + fasta + " -outfmt '6 qstart qend sstart send evalue length' > " + folder + "/" + name + ".blast"
    subprocess.run(cmd2, shell=True)
