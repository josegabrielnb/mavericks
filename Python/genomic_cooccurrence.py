#!/usr/local/bin/python3

"""This script finds co-ocurrences of two markers which are at least 4000 
and no more than 40000 base pairs apart, using tblastn result files"""

import sys, re, copy, datetime

integrase = sys.argv[-2]
polb = sys.argv[-1]

#print("Species: X")
#print("Assembly: Y")

d = {}
c = {}
delta = 40000
epsilon = 4000

with open(integrase, 'r') as f:
    llines = f.readlines()
f.closed

with open(polb, 'r') as g:
    glines = g.readlines()
g.closed

for line in llines:
    line.strip()
    if re.match("#",line) is None and re.match(r'^\s*$', line) is None:
        hits = line.split()
        region, a, b = hits[1], int(hits[8]), int(hits[9])
        if region not in d:
            d[region] = [a, b]
        else:
            d[region].append(a)
            d[region].append(b)
            
for line in glines:
    line.strip()
    if re.match("#",line) is None and re.match(r'^\s*$', line) is None:
        hits = line.split()
        region, a, b = hits[1], int(hits[8]), int(hits[9])
        if region in d:
            c[region] = d[region]
            c[region].append(a)
            c[region].append(b)
            
for key in c:
    c[key].sort()

coordinates = copy.deepcopy(c)

for key in c:
    for i in range(1, len(c[key])-1):
        step1 = c[key][i]-c[key][i-1]
        step2 = c[key][i+1]-c[key][i]
        if step1 < delta and step2 < delta:
            coordinates[key][i] = 0

for key in coordinates:
    for i in coordinates[key]:
        if i == 0:
            coordinates[key] = list(filter((0).__ne__, coordinates[key]))

candidates = copy.deepcopy(coordinates)

for key in coordinates:
    for i in range(1, len(coordinates[key])):
        step1 = coordinates[key][i]-coordinates[key][i-1]
        if step1 < epsilon:
            candidates[key][i] = 0
            candidates[key][i-1] = 0

for key in candidates:
    for i in candidates[key]:
        if i == 0:
            candidates[key] = list(filter((0).__ne__, candidates[key]))

my_result = copy.deepcopy(candidates)

for key in candidates:
    if len(candidates[key]) == 0:
        del my_result[key]

count = 0
my_sum = 0

for key in my_result:
    count += 1

for i in my_result:
    for j in my_result[i]:
        my_sum += 1

#my_results = open(species + '.out3', 'w')

now = f"{datetime.datetime.now():%d-%h-%Y %H:%M}"
print("\n",now,"\n")
#my_results.writelines(["\n", now, "\n\n"])
print("Species:","Genome assembly:", "\n")
#my_results.writelines(["%s %4s %20s %4s" % ("Species:","Genome assembly:"), "\n\n"])
print(str(delta))
#my_results.writelines(["%s %s" % ("\u0394 =", str(delta)), "\n"])
print(str(epsilon), "\n")
#my_results.writelines(["%s %5s" % ("\u03B5 =", str(epsilon)), "\n\n"])
print("# Contigs with integrase/PolB:", str(count))
#my_results.writelines(["%s %s" % ("# Contigs with integrase/PolB:", str(count)), "\n"])
print("# Candidate regions:", str(int(my_sum/2)), "\n")
#my_results.writelines(["%s %12s" % ("# Candidate regions:", str(int(my_sum/2))), "\n\n"])

for key in my_result:
    print(key, "\t", str(my_result[key]))
    #my_results.writelines([key, "\t", str(my_result[key]), "\n"])

#my_results.close()


