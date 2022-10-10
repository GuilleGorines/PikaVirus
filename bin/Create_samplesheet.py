#!/usr/bin/env python

import os
import argparse
import sys
import re

parser = argparse.ArgumentParser(description="Generates the sample csv needed for nf-core pikavirus to function.")

parser.add_argument("-directory", default="all", dest="directory", required=True, help="Directory where samples are placed, and where the samplesheet will be generated.")
args = parser.parse_args()

def find_longest_match(list_string1,list_string2):
    
    punctuation_dict = {}
    
    for string1 in list_string1:
        
        punctuation_dict[string1] = []
        dissected_1 = [letter for letter in string1]
        
        for string2 in list_string2:
            
            dissected_2 = [letter for letter in string2]
            samplename = ""
            
            for position in range(0, min(len(dissected_1),len(dissected_2))):
            
                if dissected_1[position] == dissected_2[position]:
                    samplename += dissected_1[position]
                else:
                    punctuation_dict[string1].append([string2,samplename])
                    break
    
    return punctuation_dict


def find_best_match(punctuation_dict):
    
    final_groups = []
    for key,value in punctuation_dict.items():

        secondstring = ""
        samplename = ""
        punctuation = 0

        for second,shared in value:
            if len(shared) > punctuation:
                secondstring = second
                samplename = shared
                punctuation = len(shared)
                
            else:
                pass

        # remove last R in the end so it doesnt look bad
        samplename = re.sub("_R$","",samplename)    

        final_groups.append([key,secondstring,samplename])
        
    return final_groups

def create(final_groups, path):

    truepath = os.path.realpath(path)
    
    with open(f"{truepath}/samplesheet.csv","w") as outfile:
            
        outfile.write("sample,fastq_1,fastq_2\n")
    
        for item in final_groups:
            if len(item) == 2:
                file = item[0]
                samplename = item[1]
                outfile.write(f"{samplename},{truepath}/{file},\n")
            
            elif len(item) == 3:
                file1 = item[0]
                file2 = item[1]
                samplename = item[2]
                outfile.write(f"{samplename},{truepath}/{file1},{truepath}/{file2}\n")
   
    return     
    

single_end_files = []
R1_files = []
R2_files = []

for item in os.listdir(args.directory):
    if item.endswith(".fastq") or item.endswith(".fastq.gz") or item.endswith(".fq") or item.endswith(".fgz"):
        if "_R1" in item or "_1." in item:
            R1_files.append(item)
        elif "_R2" in item or "_2." in item:
            R2_files.append(item)
        else:
            single_end_files.append(item)

# exit if different length in R1 and R2
if len(R1_files) != len(R2_files):
    print(f"There is a different number of R1 ({len(R1_files)}) and R2 ({len(R1_files)}) files.")
    sys.exit()

print(f"{len(R1_files)} paired sequences")
print(f"{len(single_end_files)} single sequences")

punctuation_dict = find_longest_match(R1_files,R2_files)
paired_files_with_samplename = find_best_match(punctuation_dict)

single_files_with_samplename = []

for single_end_file in single_end_files:
    samplename = single_end_file.split(".")[0]
    single_files_with_samplename.append([single_end_file,samplename])

final_groups = []

if len(single_files_with_samplename) != 0:
    final_groups.extend(single_files_with_samplename)

if len(paired_files_with_samplename) != 0:
    final_groups.extend(paired_files_with_samplename)

create(final_groups, args.directory)