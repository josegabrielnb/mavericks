#!/usr/local/bin/python3

"""Find the fasta sequences for candidates that were not recovered and search hits against a general
improved HMM library"""

import sys, os, subprocess, re

with open('teleostei_species.txt', 'r') as f:
    spp = f.readlines()
f.closed

for sp in spp:

    species = sp.strip()

    print("Currently analysing {}".format(species))

    os.chdir(species)

    files = os.listdir('.')

    hmms = os.listdir('/Users/user/Desktop/Genome_screening/Maverick_annotation/hmm_library2.0/')

    if '.DS_Store' in files:
        files.remove('.DS_Store')

    mavericks = []

    no_hits = re.compile(r".*No hits.*")

    for file in files:
        if re.match(r'.*results.*',file):

            with open(file,'r') as f:
                results = f.readlines()
            f.closed

            for line in results:
                my_maverick = line.split()[9]
                mavericks.append(my_maverick)

    for core in files:
        if re.match(r'Core_.*',core):

            os.chdir(core)

            fastas = os.listdir('.')

            for fasta in fastas:
                if fasta.endswith('.fa'):

                    if fasta not in mavericks:

                        dir_name = 'Recover_hmmlib2_' + core

                        recovered_files = os.listdir('..')

                        if dir_name not in recovered_files:
                            os.mkdir('../' + dir_name)

                        os.mkdir('../' + dir_name + '/' + fasta.strip())
                        
                        #ßprint(fasta)
                        #DO HMM search!

                        for hmm in hmms:
                            cmd1 = 'nhmmer --watson ~/Desktop/Genome_screening/Maverick_annotation/hmm_library2.0/' + hmm.strip() + ' ' + fasta.strip()
                            result1 = subprocess.getoutput(cmd1)
                            cmd2 = 'nhmmer --crick ~/Desktop/Genome_screening/Maverick_annotation/hmm_library2.0/' + hmm.strip() + ' ' + fasta.strip()
                            result2 = subprocess.getoutput(cmd2)

                            #print(result1)
                            #ßprint(result2)

                            if re.search(no_hits,result1) or re.search(no_hits,result2):
                                continue

                            else:
                                prefix = re.sub(r'.fa','_', fasta.strip())
                                print("Two matches found for {} in {}".format(fasta,core))

                            if hmm.endswith('5prime200.hmm'):
                            
                                watson1 = re.sub(r'5prime200.hmm','watson.txt', hmm)
                                watson = prefix + watson1

                                crick1 = re.sub(r'5prime200.hmm','crick.txt', hmm)
                                crick = prefix + crick1

                            elif hmm.endswith('5prime200_AGTAGT.hmm'):

                                watson1 = re.sub(r'5prime200_AGTAGT.hmm','AGTAGT_watson.txt', hmm)
                                watson = prefix + watson1

                                crick1 = re.sub(r'5prime200_AGTAGT.hmm','AGTAGT_crick.txt', hmm)
                                crick = prefix + crick1

                            elif hmm.endswith('5prime200_AGAGAG.hmm'):

                                watson1 = re.sub(r'5prime200_AGAGAG.hmm','AGAGAG_watson.txt', hmm)
                                watson = prefix + watson1

                                crick1 = re.sub(r'5prime200_AGAGAG.hmm','AGAGAG_crick.txt', hmm)
                                crick = prefix + crick1


                            print("Hits found between {} and {}".format(fasta.strip(),hmm.strip()))
                            cmd3 = 'nhmmer --watson --tblout ' + '../' + dir_name + '/' + fasta.strip() + '/' + watson + ' ~/Desktop/Genome_screening/Maverick_annotation/hmm_library2.0/' + hmm.strip() + ' ' + fasta.strip()
                            subprocess.run(cmd3,shell=True)
                            cmd4 = 'nhmmer --crick --tblout ' + '../' + dir_name + '/' + fasta.strip() + '/' + crick + ' ~/Desktop/Genome_screening/Maverick_annotation/hmm_library2.0/' + hmm.strip() + ' ' + fasta.strip()
                            subprocess.run(cmd4,shell=True)

            os.chdir('..')

    print("\n")
    print(mavericks)

    os.chdir('..')

