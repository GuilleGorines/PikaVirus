#!/usr/bin/env python

'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 1.1

CREATED: Exact date unknown

REVISED: 26-7-2021

DESCRIPTION: 
    Checks MASH result, extracts the name of the significative references 
    found in the given directory, and creates a symlink for Nextflow 
    to pick for future references.

INPUT (by order):
    1. MASH file result (tsv format, from mash dist)
    2. path to directory containing the reference data


USAGE:
    extract_significative_references.py MASHfile directory

REQUIREMENTS:
    -Python >= 3.6

DISCLAIMER: this script has exclusively been developed for nf-core pikavirus,
and therefore is not guaranteed to function properly in other settings. 
Despite this, feel free to use it at will.

TO DO: 

================================================================
END_OF_HEADER
================================================================
'''

import sys
import os

mashresult = sys.argv[1]
refdir = sys.argv[2]
ref_sheet = sys.argv[3]
realpath = os.path.realpath(refdir)

# 0 Identity
# 1 Shared_hashes
# 2 Median_multiplicity
# 3 P-value
# 4 Query_id
# 5 Query_comment


with open(mashresult) as infile:
    infile = infile.readlines()

    # remove header 
    infile = [line.split("\t") for line in infile if not line.startswith("#")]

    chosen = []
    # Criteria:
    # Identity of over 0.9
    # P-val over 0.05
    # More than half shared hashes
    
    for line in infile:

        if float(line[0]) > 0.9 and float(line[3]) < 0.05:

            numerator, denominator = line[1].split("/")
            shared = int(numerator)/int(denominator)

            if shared > 0.1:
                chosen.append(line[4].split("/")[-1])

# Reference name, not only the file

file_headers = ["filename","file_name","file-name","file"]
species_name_headers = ["scientific_name","organism_name","organism","species_name","species"]
subspecies_name_headers = ["intraespecific_name","subspecies_name","strain","subspecies"]

# Read sample sheet
with open(ref_sheet) as reference_data:    
    reference_data = reference_data.readlines()
    headers = [line.split("\t") for line in reference_data if line.startswith("#")]
    reference_data = [line.split("\t") for line in reference_data if not line.startswith("#")]

# Identify which columns are here
for single_header in headers:
    for item in single_header:
        if item.lower() in file_headers:
            file_column_index = single_header.index(item)

        elif item.lower() in species_name_headers:  
            species_column_index = single_header.index(item)

        elif item.lower() in subspecies_name_headers:
            subspecies_column_index = single_header.index(item)


reference_data = [line for line in reference_data if line[file_column_index] in chosen]

species_dict = {}

# put the different results in a dict, shared if no strain

for item in reference_data:
    if item[subspecies_column_index] == "":
        entry = item[species_column_index]
    else:
        entry = f"{item[species_column_index]} {item[subspecies_column_index]}"

    size = os.path.getsize(f"{realpath}/{item[file_column_index]}")
    
    # get only the biggest file as representative
    if entry not in species_dict.keys():
        species_dict[entry] = [item[file_column_index], size]
    else:
        if size > species_dict[entry][1]:
            species_dict[entry] = [item[file_column_index], size]

print(species_dict)

chosen = [item[0] for item in species_dict.values()]  

# files
reference_dict = {item:[f"{realpath}/{item}",f"Final_fnas/{item}"] for item in os.listdir(refdir) if item in chosen} 

os.mkdir(f"Final_fnas", 0o777)

for assembly in chosen:
    os.symlink(f"{realpath}/{assembly}",f"Final_fnas/{assembly}")

# If no coincidences:
if not os.listdir("Final_fnas"):
    with open("not_found.tsv","w") as outfile:
        outfile.write("NO ORGANISMS FOUND")
 
