#!/usr/local/bin/python3

"""Identify terminal inverted repeats in a database of nucleotide sequence
which were blasted against themselves"""

import sys, os, re

folder = sys.argv[-1]

files = os.listdir(folder)

my_files = [my_file for my_file in files if bool(re.search(".blast", my_file))]

for my_file in my_files:

    name = my_file.replace(".blast","")

    path = folder + "/" + my_file

    with open(path, 'r') as f:
        lines = f.readlines()
    f.closed

    for line in lines:
        line = line.split()
        a, b, c, d = int(line[0]), int(line[1]), int(line[2]), int(line[3])
        tir_length = int(line[5])
        maverick_length = abs(d - a)
        if (c > d) and (tir_length > 80) and (4000 <= maverick_length <= 40000):
            print(name, a, b, c, d, tir_length, line[4], maverick_length)