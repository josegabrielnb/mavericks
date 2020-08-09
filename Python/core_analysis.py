#!/usr/local/bin/python3

"""This script analyses the blast hit table of the core Maverick genes to score presence/absence"""

import sys, os, subprocess

species = sys.argv[-1]

with open(species, 'r') as f:
    spp = f.readlines()
f.closed

for sp in spp:
    folder = sp.rstrip()
    for blast in os.listdir(folder):
        if blast.endswith(".blast"):
            print("Analysing core genes for {} ({})...".format(blast,folder))
            name = blast.replace(".blast","")
            table = {}
            count = 0

            path1 = folder + "/" + blast

            with open(path1,'r') as g:
                lines = g.readlines()
            g.close()

            for line in lines:
                datum = line.split()
                marker = datum[0]
                evalue = float(datum[2])
                if marker not in table.keys():
                    table[marker] = evalue
                    if evalue < 1e-5:
                        count += 1

            path2 = folder + "/" + name + ".core"

            outfile = open(path2, 'w')
            title = "BEST HITS BELOW THRESHOLD (e-value < 1e-5): " + str(count) + "\n\n"
            outfile.write(title)
            for key, value in table.items():
                    out = str(key) + "\t" + str(value) + "\n"
                    outfile.write(out)
            outfile.close()
            print("Done!")
