#!/usr/local/bin/python3

"""This script finds target site duplications in nucleotide sequences"""

import sys, re

species = sys.argv[-1]

with open(species, 'r') as f:
    lines = f.readlines()
f.closed

seq = ''
motifs = []
n = 0
start_search1 = 14000
end_search1 = 14712
start_search2 = 30654
end_search2 = 31800

for line in lines:
    if re.match(">",line) is None:
        line = line.strip()
        seq = seq + line

# function by Tommy Tang
def ReverseComplement1(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

for i in range(start_search1,end_search1-6):
    tsd1 = seq[i:i+6]
    end5 = seq[i+6:i+12]
    for j in range(start_search2,end_search2-6):
        tsd2 = seq[j:j+6]
        end3 = seq[j-6:j]
        revcomp3 = ReverseComplement1(end3)
        motif_score = 0
        tsd_score = 0
        for k in range(0,6):
            if end5[k] == revcomp3[k]:
                motif_score += 1
            if tsd1[k] == tsd2[k]:
                tsd_score += 1
        score = motif_score + tsd_score
        if score >= 10:
                print("{} {} {} {} {} {}-{} Total Score {}".format(tsd1, tsd2 ,end5, end3, revcomp3, i+1 ,j+6, score))
                n += 1

print(n)