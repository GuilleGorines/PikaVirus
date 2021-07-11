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

REVISED: 26-5-2021

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

with open(mashresult) as infile:
    infile = infile.readlines()

    # remove header 
    infile = [line.split() for line in infile if not line.startswith("#")]

    # get name of file if p-val < 0.05
    infile = [line[0].split("/")[-1] for line in infile if float(line[3]) < 0.05]
    
# files
reference_dict = {item.split()[0]:[f"{realpath}/{item}",f"Final_fnas/{item}"] for item in os.listdir(refdir)}

os.mkdir(f"Final_fnas", 0o777)

for filename in infile:
    if filename in reference_dict.keys():
        os.symlink(reference_dict[filename][0],reference_dict[filename][1])

# If no coincidences:
if not os.listdir("Final_fnas"):
    with open("not_found.txt","w") as outfile:
        outfile.write("#Header")
        outfile.write("NO ORGANISMS FOUND")
 
