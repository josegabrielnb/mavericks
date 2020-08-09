#!/usr/local/bin/python3

"""nhmmer against a HMM library for all species iteration 2""""

import sys, os, subprocess, re

hmms = os.listdir('/Users/user/Desktop/Genome_screening/Maverick_annotation/hmm_library')
cores = os.listdir('.')
if '.DS_Store' in cores:
    cores.remove('.DS_Store')

no_hits = re.compile(r".*No hits.*")

for core in cores:
    if re.match(r'Core',core):
        os.chdir(core)
        files = os.listdir('.')

        for file in files:
            if file.endswith('.fa'):

                dir_name = 'Recover_hmmlib_' + core

                recovered_files = os.listdir('..')

                if dir_name not in recovered_files:
                    os.mkdir('../' + dir_name)

                os.mkdir('../' + dir_name + '/' + file.strip())

                for hmm in hmms:
                    cmd1 = 'nhmmer --watson ~/Desktop/Genome_screening/Maverick_annotation/hmm_library/' + hmm.strip() + ' ' + file.strip()
                    result1 = subprocess.getoutput(cmd1)
                    cmd2 = 'nhmmer --crick ~/Desktop/Genome_screening/Maverick_annotation/hmm_library/' + hmm.strip() + ' ' + file.strip()
                    result2 = subprocess.getoutput(cmd2)
                    if re.search(no_hits,result1) or re.search(no_hits,result2):
                        continue
                    else:

                        prefix = re.sub(r'.fa','_', file.strip())

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


                        print("Hits found between {} and {}".format(file.strip(),hmm.strip()))
                        cmd3 = 'nhmmer --watson --tblout ' + '../' + dir_name + '/' + file.strip() + '/' + watson + ' ~/Desktop/Genome_screening/Maverick_annotation/hmm_library/' + hmm.strip() + ' ' + file.strip()
                        subprocess.run(cmd3,shell=True)
                        cmd4 = 'nhmmer --crick --tblout ' + '../' + dir_name + '/' + file.strip() + '/' + crick + ' ~/Desktop/Genome_screening/Maverick_annotation/hmm_library/' + hmm.strip() + ' ' + file.strip()
                        subprocess.run(cmd4,shell=True)

        os.chdir('..')
