#!/usr/local/bin/python3

"""This script uses HMMs of the first 200-bp of the 5' inverted repeat to recover the coordinates
of AGAGAG mavericks that could not be detected with maverick_miner.py"""

import sys, os, subprocess, re

with open("perca_species.txt", 'r') as f:
            spp = f.readlines()
f.closed

for sp in spp:

    species = sp.strip()

    print("Analysing {} for additional Mavericks...".format(species))

    os.chdir(species)

    with open(species + "_results.txt", 'r') as f:
        lines = f.readlines()
    f.closed

    #cores = next(os.walk('.'))[1]

    seqs = {}
    recovery_seqs = {}

    files = os.listdir('.')

    for file in files:
        if re.match(r'Core', file):
                seqs[file] = []
                recovery_seqs[file] = []

    for line in lines:
        data = line.rstrip().split()
        seqs[data[1]].append(data[9])

    for key in seqs.keys():
        os.chdir(key)

        files = os.listdir('.')

        for file in files:
            if file.endswith('.fa') and file not in seqs[key]:
                recovery_seqs[key].append(file)

        os.chdir('..')

    print(recovery_seqs)

    for core in recovery_seqs.keys():

        dir_name = 'Recover_AGAGAG_' + core

        os.mkdir(dir_name)

        for fasta in recovery_seqs[core]:

            table_name1 = re.sub(r'.fa','_watson.txt', fasta)
            cmd1 = 'nhmmer --noali --watson --tblout ' + dir_name + '/' + table_name1 + ' ' + species + '_5prime200_AGAGAG.hmm ' + core + '/' + fasta
            subprocess.run(cmd1, shell = True)
            table_name2 = re.sub(r'.fa','_crick.txt', fasta)
            cmd2 = 'nhmmer --noali --crick --tblout ' + dir_name + '/' + table_name2 + ' ' + species + '_5prime200_AGAGAG.hmm '  + core + '/' +  fasta
            subprocess.run(cmd2, shell = True)

    #EXTRACT CANDIDATE COORDINATES AND CHECK TSDs

    for core in recovery_seqs.keys():

        dir_name = 'Recover_AGAGAG_' + core

        for fasta in recovery_seqs[core]:

            table_name1 = re.sub(r'.fa','_watson.txt', fasta)
            table_name2 = re.sub(r'.fa','_crick.txt', fasta)

            with open(dir_name + '/' + table_name1, 'r') as w:
                watson = w.readlines()
            w.closed

            with open(dir_name + '/' + table_name2, 'r') as c:
                crick = c.readlines()
            c.closed

            end5 = watson[2].split()
            end3 = crick[2].split()

            if len(end5) > 1 and len(end3) > 1:
                start = int(end5[6]) - int(end5[4]) + 1
                end = int(end3[6]) + int(end3[4]) - 1
                
                if int(end5[4]) < 10 and int(end3[4]) < 10:

                    #FIND MOTIFS OF TIRs AND TSDs

                    cmd3 = 'blastdbcmd -db ' + core + '/' + fasta + ' -entry 1 -range ' + str(start-6) + '-' + str(start+5)
                    result5 = subprocess.getoutput(cmd3)

                    cmd4 = 'blastdbcmd -db ' + core + '/' + fasta + ' -entry 1 -range ' + str(end-5) + '-' + str(end+6)
                    result3 = subprocess.getoutput(cmd4)

                    seq5 = result5.split('\n')
                    seq3 = result3.split('\n')

                    if (end-start) > 0 and abs(end-start) < 50000:
                        #print(start,end,end-start,seq5[1][0:6],seq5[1][6:12],seq3[1][0:6],seq3[1][6:12])
                        cmd5 = 'blastdbcmd -db ' + core + '/' + fasta + ' -entry 1 -range ' + str(start) + '-' + str(end)
                        sequence = subprocess.getoutput(cmd5)

                        loc1 = re.sub(r'.fa', '', fasta)
                        loc2 = re.sub(r'.*[1-9]\.', '', loc1)

                        nums = loc2.split("-")

                        if nums[0] == '':
                            start2 = 1
                        else:
                            start2 = int(nums[0])

                        name = ">" + re.sub(r'\.-*[0-9]*-[0-9]*\.fa', ':', fasta) + str(start + start2 - 1) + "-" + str(end + start2 - 1)
                        curated_sequence = sequence.replace(">1", name)

                        AGT_motif = 'AGTAGT'
                        ACT_motif = 'ACTACT'
                        AG_motif = 'AGAGAG'
                        CT_motif = 'CTCTCT'
                        AGT_score = 0
                        ACT_score = 0
                        AG_score = 0
                        CT_score = 0

                        for i in range(0,6):
                            if seq5[1][i+6] == AGT_motif[i]:
                                AGT_score += 1
                            elif seq5[1][i+6] == AG_motif[i]:
                                AG_score += 1

                        for i in range(0,6):
                            if seq3[1][i] == ACT_motif[i]:
                                ACT_score += 1
                            elif seq3[1][i] == CT_motif[i]:
                                CT_score += 1

                        AGT_total = AGT_score + ACT_score
                        AG_total = AG_score + CT_score

                        if AGT_total > AG_total:
                            with open(species + "_recovered_mavericks_AGTAGT.fa", 'a+') as o:
                                print(curated_sequence, file=o)
                            o.closed
                        else:
                            with open(species + "_recovered_mavericks_AGAGAG.fa", 'a+') as o:
                                print(curated_sequence, file=o)
                            o.closed

                        table = species + '_recovered_results.txt'

                        with open(table, 'a+') as f:
                            print(species,core,seq5[1][0:6],seq5[1][6:12],seq3[1][0:6],seq3[1][6:12],start + start2 - 1,end + start2 - 1,end-start,fasta,file=f)
                        f.closed

                        print(seq5[1][0:6],seq5[1][6:12],AGT_total,AG_total,seq3[1][0:6],seq3[1][6:12])

    os.chdir('..')
            #table_name2 = re.sub(r'.fa','_crick.txt', fasta)