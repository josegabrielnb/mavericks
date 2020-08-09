#!/usr/local/bin/python3

"""This script analyses the blast hit table of the core Maverick genes and counts best hits
below the evalue threshold"""

import sys, os, subprocess

species = sys.argv[-1]

with open(species, 'r') as f:
    spp = f.readlines()
f.closed

for sp in spp:
    folder = sp.rstrip()
    core_table = {'Core_0':0,'Core_1':0,'Core_2':0,'Core_3':0,'Core_4':0,'Core_5':0,'Core_6':0}
    core_seqs = {'Core_0':[],'Core_1':[],'Core_2':[],'Core_3':[],'Core_4':[],'Core_5':[],'Core_6':[]}
    for core in os.listdir(folder):
        if core.endswith(".core"):

            name = core.replace(".core",".fa")
            path1 = folder + "/" + core

            with open(path1,'r') as g:
                lines = g.readlines()
            g.close()

            core_genes = lines[0].rstrip().split()

            tag = "Core_" + core_genes[7]
            core_table[tag] += 1
            core_seqs[tag].append(name)

    print(folder)
    for k,v in core_table.items():
        print("{} {}".format(k,v))
    for k,v in core_seqs.items():
        print("{} {}".format(k,v))

