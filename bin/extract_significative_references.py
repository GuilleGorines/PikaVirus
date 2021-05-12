#!/usr/bin/env python
#
# USAGE:
# Checks MASH result, extracts the name of the significative references found in the given directory,
# and creates a symlink for nextflow to pick for future references.
# 
# INPUTS:
# - 1: MASH file result (tsv format, obtained from mash dist)
# - 2: path to directory containing the reference data
#
#
# DISCLAIMER: this script has exclusively been developed for the correct functioning of nf-core pikavirus,
# and therefore is not guaranteed to function properly in other settings. Despite this, feel free to use it
# at will.

import sys
import os

mashresult = sys.argv[1]
refdir = sys.argv[2]
realpath = os.path.realpath(refdir)

with open(mashresult) as infile:
    infile = infile.readlines()

    # remove header 
    infile = [line.split() for line in infile if not line.startswith("#")]

    # get name, remove extension of file if p-val < 0.05
    infile = [line[0].split()[0] for line in infile if float(line[3]) < 0.05]
    
# files without extension 
reference_dict = {item.split()[0]:[f"{realpath}/{item}",f"Final_fnas/{item}"] for item in os.listdir(refdir)}

os.mkdir(f"Final_fnas", 0o777)

for filename in infile:
    if filename in reference_dict.keys():
        os.symlink(reference_dict[filename][0],reference_dict[filename][1])


