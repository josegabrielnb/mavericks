#!/usr/local/bin/python3

"""Make phase 3 hmms for remaining fish species (single Maverick type)"""

import sys, os, subprocess, re

with open('teleostei_species_phase3.txt','r') as f:
    spp = f.readlines()
f.closed

for sp in spp:
    species = sp.strip()

    os.chdir(species)

    files = os.listdir('.')

    a = ''
    b = ''
    c = ''

    for file in files:
        if file == species + '_mavericks.fa':
            a = file.strip()
    
    for file in files:
          if file == species + '_recovered_mavericks_AGAGAG.fa':
            b = file.strip()

    for file in files:
          if file == species + '_recovered_mavericks_AGTAGT.fa':
            c = file.strip()

    print("Working on {}".format(species))

    cmd1 = 'cat ' + a + ' ' + b + ' ' + c + ' > ' + species + '_phase12_mavericks.fa'
    subprocess.run(cmd1,shell=True)
    print("Concatenation successful!")
    cmd2 = 'mafft ' + species + '_phase12_mavericks.fa > ' + species + '_phase12_mavericks.mafft.aln'
    subprocess.run(cmd2,shell=True)
    print("mafft alignment successful!")
    cmd3 = 'extractalign -sequence ' + species + '_phase12_mavericks.mafft.aln' + ' -regions 1-200 -outseq ' + species + '_phase12_5prime200.aln'
    subprocess.run(cmd3,shell=True)
    print("5' prime extraction successful!")
    cmd4 = 'hmmbuild ../../hmm_library2.0/' + species + '_phase12_5prime200.hmm' + ' ' + species + '_phase12_5prime200.aln'
    subprocess.run(cmd4,shell=True)
    print("hmmbuild successful!\n")

    os.chdir('..')