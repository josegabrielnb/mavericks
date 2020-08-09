#!/usr/local/bin/python3

""""Analyse data for files produced by hmm_species_search.py"""

import sys, os, subprocess, re

with open("teleostei_species_iteration2.txt", 'r') as f:
            spp = f.readlines()
f.closed

for sp in spp:

    species = sp.strip()

    os.chdir(species)

    directories = os.listdir('.')

    for directory in directories:

        if re.match(r'Recover_hmmlib*',directory):
            #print(directory)
            os.chdir(directory)
            fastas = os.listdir('.')

            for fasta in fastas:

                if '.DS_Store' not in fasta:
                    
                    #print("\n" + fasta)

                    os.chdir(fasta.strip())
                    files = os.listdir('.')

                    hmm_5end = {}
                    hmm_3end = {}

                    for file in files:
                        
                        if re.match(r'.*watson.txt',file):

                            with open(file,'r') as w:
                                watson = w.readlines()
                            w.closed
                            
                            for line in watson:
                                if '#' not in line:
                                    data = line.split()

                                    if not hmm_5end:
                                        hmm_5end[fasta] = [data[2],data[4],data[6],data[12]]
                                    elif float(hmm_5end[fasta][3]) > float(data[12]):
                                        hmm_5end[fasta] = [data[2],data[4],data[6],data[12]]

                        elif re.match(r'.*crick.txt',file):

                            with open(file,'r') as c:
                                crick = c.readlines()
                            c.closed

                            for line in crick:
                                if '#' not in line:
                                    data = line.split()

                                    if not hmm_3end:
                                        hmm_3end[fasta] = [data[2],data[4],data[6],data[12]]
                                    elif float(hmm_3end[fasta][3]) > float(data[12]):
                                        hmm_3end[fasta] = [data[2],data[4],data[6],data[12]]

                    #print(hmm_5end)
                    #print(hmm_3end)

                    if hmm_5end and hmm_3end:

                        if int(hmm_5end[fasta][1]) < 10 and int(hmm_3end[fasta][1]) < 10:

                            folder = re.sub(r'Recover_hmmlib_','',directory)
                            end5 = int(hmm_5end[fasta][2]) - int(hmm_5end[fasta][1]) + 1
                            end3 = int(hmm_3end[fasta][2]) + int(hmm_3end[fasta][1]) - 1

                            if end5 < end3:

                                cmd1 = 'blastdbcmd -db ' + '../../' + folder + '/' + fasta.strip() + ' -entry 1 -range ' + str(end5) + '-' + str(end3)
                                sequence = subprocess.getoutput(cmd1)
                                cmd2 = 'blastdbcmd -db ' + '../../' + folder + '/' + fasta.strip() + ' -entry 1 -range ' + str(end5-6) + '-' + str(end5+5)
                                cmd3 = 'blastdbcmd -db ' + '../../' + folder + '/' + fasta.strip() + ' -entry 1 -range ' + str(end3-5) + '-' + str(end3+6)
                                seq5 = subprocess.getoutput(cmd2)
                                seq3 = subprocess.getoutput(cmd3)

                                seq5 = seq5.split("\n")
                                seq3 = seq3.split("\n")

                                loc1 = re.sub(r'.fa', '', fasta)
                                loc2 = re.sub(r'.*[1-9]\.', '', loc1)

                                nums = loc2.split("-")

                                if nums[0] == '':
                                    start2 = 1
                                else:
                                    start2 = int(nums[0])

                                name = ">" + re.sub(r'\.-*[0-9]*-[0-9]*\.fa', ':', fasta) + str(end5 + start2 - 1) + "-" + str(end3 + start2 - 1)
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
                                    with open('../../' + species + "_recovered_mavericks_AGTAGT.fa", 'a+') as o:
                                        print(curated_sequence, file=o)
                                    o.closed
                                else:
                                    with open('../../' + species + "_recovered_mavericks_AGAGAG.fa", 'a+') as o:
                                        print(curated_sequence, file=o)
                                    o.closed

                                table = '../../' + species + '_recovered_results.txt'

                                with open(table, 'a+') as f:
                                    print(species,folder,seq5[1][0:6],seq5[1][6:12],seq3[1][0:6],seq3[1][6:12],str(end3-end5),str(end5 + start2 - 1),str(end3 + start2 - 1),fasta,file=f)
                                f.closed

                                print(curated_sequence)
                                print(seq5[1][0:6],seq5[1][6:12],seq3[1][0:6],seq3[1][6:12])

                    os.chdir('..')
        
            os.chdir('..')

    os.chdir('..')