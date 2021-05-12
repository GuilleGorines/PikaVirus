#!/usr/bin/env python

# USAGE:
# Checks kraken2 report and extracts those files in the reference that respond to the species (S) taxids.
# 
# INPUTS:
# - 1: kraken report
# - 2: tsv file containing the filename and the asssociated taxid
# - 3: directory containing the reference mentioned in input 2

#    NOTE: this file must refer to files given to the corresponding reference directory, and
#    can be automatically generated in the scripts given by nf-core/pikavirus. Otherwise, the 
#    file schema MUST be:
#
#    *headers not necessary*
#    col0: Assembly_accession
#    col1: Species_taxid
#    col2: Subspecies_taxid (= to Species_taxid if not a subspecies)
#    col3: Scientific_name
#    col4: Intraespecific_name
#    col5: Download_URL
#    col6: File_name
#
#   This schema is temporal, and will probably be updated to be less-restrictive requirement in 
#   future releases. Don't worry, we are actively working on it!
#
# DISCLAIMER: this script has exclusively been developed for the correct functioning of nf-core pikavirus,
# and therefore is not guaranteed to function properly in other settings. Despite this, feel free to use it
# at will.

import sys
import os
import shutil

krakenrep = sys.argv[1]
reference_naming = sys.argv[2]
reference_directory = sys.argv[3]

realpath = os.path.realpath(reference_directory)
# Report schema:
#   3: rank code (Unclass, Kingdom...)
#   4: taxID
#   5: scientific_name (indented)

# Look for "S" or "S1" in the report and extract Taxid
with open(krakenrep) as krakenreport:
    krakenreport = [line.split("\t") for line in krakenreport.readlines()]
    idlist = []
    for line in krakenreport:
        if line[3] == "S" or line[3] == "S1":
            idlist.append(line[4])

    if len(idlist) == 0:
        print(f"No species were identified in the kraken report.")
        sys.exit(2)    

# Look for the found taxids in the reference file:
with open(reference_naming) as refids:
    refids = refids.readlines()
    refids = [line.strip("\n").split("\t") for line in refids]
    refids = [line for line in refids if line[1] in idlist or line[2] in idlist]

# Extract filenames
filelist = [line[6] for line in refids]

# Look for present filenames path
filedict = {filename:[f"{realpath}/{filename}",f"Chosen_fnas/{filename}"] for filename in os.listdir(reference_directory)}


os.mkdir(f"Chosen_fnas", 0o777)

for filename in filedict.keys():
    if filename in filelist:
        os.symlink(filedict[filename][0], filedict[filename][1])