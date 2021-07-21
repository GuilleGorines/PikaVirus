#!/usr/bin/env python

'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 1.0

CREATED: Exact date unknown

REVISED: 20-7-2021

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

            if shared > 0.5:
                chosen.append(line[4].split("/")[-1])
    
# files
reference_dict = {item:[f"{realpath}/{item}",f"Final_fnas/{item}"] for item in os.listdir(refdir) if item in chosen} 

os.mkdir(f"Final_fnas", 0o777)

for assembly in chosen:
    os.symlink(f"{realpath}/{assembly}",f"Final_fnas/{assembly}")

# If no coincidences:
if not os.listdir("Final_fnas"):
    with open("not_found.tsv","w") as outfile:
        outfile.write("#Header\n")
        outfile.write("NO ORGANISMS FOUND")
 
