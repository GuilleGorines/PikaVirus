#!/usr/bin/env python

'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 1.0 

CREATED: Exact date unknown (about early 2021)

REVISED: 26-5-2021

DESCRIPTION: 
    checks Kraken 2 report, and extracts those files in the reference directory 
    that respond to the species(S) or subspecies(S1), according to the ref tsv
    given

INPUTS (by order):
    1. Kraken 2 report
    2. TSV file containing the filename and the asssociated taxid
        Please check format in the pikavirus manual if using independently
    3. Directory containing the reference mentioned in input 2

USAGE:
    reference_choosing.py kraken2_report

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
    headers = [line.strip("\n").split("\t") for line in refids if line.startswith("#")]
    refids = [line.strip("\n").split("\t") for line in refids if not line.startswith("#")]
    
file_headers = ["filename","file_name","file-name","file"]
species_taxid_headers = ["species_taxid","organism_taxid"]
subspecies_taxid_headers = ["subspecies_taxid","strain_taxid"]

# find headers corresponding to filename, speciesID and subspeciesID
for single_header in headers:
    for item in single_header:
        if item.lower() in file_headers:
            file_column_index = single_header.index(item)

        elif item.lower() in species_taxid_headers:
            species_column_index = single_header.index(item)

        elif item.lower() in subspecies_taxid_headers:
            subspecies_column_index = single_header.index(item)

# Exit with error status if one of the required groups is not identified
if not file_column_index or not species_column_index or not subspecies_column_index:
    if not file_column_index:
        print(f"No headers indicating \"File name\" were found.")
    if not species_column_index:
        print(f"No headers indicating \"Species taxid\" were found.")
    if not subspecies_column_index:
        print(f"No headers indicating \"Subspecies taxid\" were found.")

    print(f"Please consult the reference sheet format, sorry for the inconvenience!")
    sys.exit(1)



# Extract filenames
filelist = [line[file_column_index] for line in refids if line[species_column_index] in idlist or line[subspecies_column_index] in idlist]

# Remove extensions
file_extensions = [".fna",".gz"]

filelist_noext = []

for item in filelist:
    item_noext = item
    for extension in file_extensions:
        item_noext = item_noext.replace(extension,"")
    filelist_noext.append(item_noext)

os.mkdir(f"Chosen_fnas", 0o777)

for filename in os.listdir(reference_directory):
    
    filename_noext = filename
    for extension in file_extensions:
        filename_noext = filename_noext.replace(extension,"")

    if filename_noext in filelist_noext:
        
        origin = f"{realpath}/{filename}"
        destiny = f"Chosen_fnas/{filename}"

        os.symlink(origin, destiny)