#!/usr/local/bin/python3

"""This script finds the genomic coordinates of complete Maverick elements and the conserved
terminal motifs if they are associated with target site duplications"""

import sys, os, subprocess, re

with open("teleostei_species.txt", 'r') as f:
            spp = f.readlines()
f.closed

for sp in spp:

    print("\nAnalysing {} for complete Mavericks...\n".format(sp.strip()))

    os.chdir(sp.strip())

    folders = os.listdir('.')

    if ".DS_Store" in folders:
        folders.remove('.DS_Store')

    for folder in folders:

        os.chdir(folder.strip())

        files = os.listdir('.')

        my_mavericks = {}

        for file in files:
            if file.endswith(".inverted"):
                with open(file, 'r') as f:
                    lines = f.readlines()
                f.closed
                n = len(lines)
                for n in range(0,len(lines)):
                    lines[n].rstrip()
                    if 'AGTAGT' or 'AGAGAG' in lines[n]:
                        data1 = lines[n].split()
                        data2 = lines[n-1].split()
                        end5 = min(int(data1[0]),int(data1[1]),int(data1[2]),int(data1[3]),int(data2[0]),int(data2[1]),int(data2[2]),int(data2[3]))
                        end3 = max(int(data1[0]),int(data1[1]),int(data1[2]),int(data1[3]),int(data2[0]),int(data2[1]),int(data2[2]),int(data2[3]))

                        fasta = file.replace('inverted','fa')

                        cmd3 = 'blastdbcmd -db ' + fasta + ' -entry 1 -range ' + str(end5-20) + '-' + str(end3+20)
                        #Right directory!
                        result = subprocess.getoutput(cmd3)

                        raw_seq = result.split("\n")
                        my_seq = ''
                        for line in raw_seq:
                            if '>' not in line:
                                my_seq += line

                        for i in range(0,50):
                            tsd1 = my_seq[i:i+6]
                            terminus_5 = my_seq[i+6:i+12]
                            for j in range(len(my_seq)-50,len(my_seq)):
                                tsd2 = my_seq[j:j+6]
                                terminus_3 = my_seq[j-6:j]

                                if tsd1 == tsd2 and (terminus_5 == 'AGTAGT' or terminus_5 == 'AGAGAG') and (terminus_3 == 'ACTACT' or terminus_3 == 'CTCTCT'):
                                    features = [end5-20+i+6, end5-20+j-1, tsd1, terminus_5, terminus_3, tsd2, j-i-6]
                                    
                                    if fasta not in my_mavericks.keys():
                                        my_mavericks[fasta] = features
                                        cmd4 = 'blastdbcmd -db ' + fasta + ' -entry 1 -range ' + str(features[0]) + '-' + str(features[1])
                                        sequence = subprocess.getoutput(cmd4)

                                        loc1 = re.sub(r'.fa', '', fasta)
                                        loc2 = re.sub(r'.*[1-9]\.', '', loc1)

                                        nums = loc2.split("-")

                                        if nums[0] == '':
                                            start = 1
                                        else:
                                            start = int(nums[0])

                                        name = ">" + re.sub(r'\.-*[0-9]*-[0-9]*\.fa', ':', fasta) + str(start + features[0] - 1) + "-" + str(start + features[1] - 1)
                                        curated_sequence = sequence.replace(">1", name)

                                        my_mavericks[fasta].extend([start + features[0] - 1, start + features[1] - 1])

                                        print("Complete Maverick found in {}".format(fasta))

                                        out1 = '../' + sp.strip() + '_complete_' + folder + '_AGAGAG' + '.fa'
                                        out2 = '../' + sp.strip() + '_complete_' + folder + '_AGTAGT' + '.fa'
                                        
                                        if terminus_5 == 'AGAGAG':
                                            with open(out1, 'a+') as f:
                                                print(curated_sequence, file = f)
                                            f.closed
                                        elif terminus_5 == 'AGTAGT':
                                            with open(out2, 'a+') as f:
                                                print(curated_sequence, file = f)
                                            f.closed

        table = '../' + sp.strip() + '_results.txt'

        with open(table, 'a+') as f:
            for maverick in my_mavericks.keys():
                print(sp.strip(), folder, *my_mavericks[maverick][2:], maverick, sep = " ", file = f)
        f.closed

        os.chdir('..')

    print("Done!\n")
        
    os.chdir('..')
