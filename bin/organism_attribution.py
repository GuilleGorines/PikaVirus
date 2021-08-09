#!/usr/bin/env python

'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 1.0

CREATED: 28-7-2021

REVISED: 28-7-2021

DESCRIPTION:
    Organizes the consensus sequences in directories, so each directory contains files for the same species (strains included)

INPUT (by order):
    -Sample name

USAGE:

REQUIREMENTS:
    -Python >= 3.6

DISCLAIMER:

TO DO:

================================================================
END_OF_HEADER
================================================================
'''

import os
import sys

samplename  = sys.argv[1]
data_sheet = sys.argv[2]

consensus_sequences = sys.argv[3:]

sequence_references = []

for item in consensus_sequences:
    
    # retrieve the reference name used 
    name = item.split("_organism_")[-1].replace(".mpileup","").replace("[","").replace("]","").replace("_consensus.fa","")
    sequence_references.append([item,name])

# Parse reference data
with open(data_sheet) as reference_data:
    reference_data = reference_data.readlines()
    headers = [line.split("\t") for line in reference_data if line.startswith("#")]
    reference_data = [line.split("\t") for line in reference_data if not line.startswith("#")]

# identify headers index
file_headers = ["filename","file_name","file-name","file"]
species_name_headers = ["scientific_name","organism_name","organism","species_name","species"]

for single_header in headers:
    for item in single_header:
        if item.lower() in file_headers:
            file_column_index = single_header.index(item)

        elif item.lower() in species_name_headers:  
            species_column_index = single_header.index(item)

# extract filename and species
filenames_in_datasheet = [[item[file_column_index], item[species_column_index]] for item in reference_data]

consensus_dict = {}

# consensus_sequence = name of the fasta containing the consensus sequence
# reference_name_consensus = name of the reference file used to extract the consensus sequence 

# reference_name_datasheet = name of (each) reference file in the datasheet
# species = name of the species 

for consensus_sequence, reference_name_consensus in sequence_references:
    
    # if filename == name in ref sheet (reference_name_datasheet), add to a dict with spp as key
    for reference_name_datasheet, species in filenames_in_datasheet:

        if reference_name_datasheet == reference_name_consensus:

            if species in consensus_dict.keys():
                consensus_dict[species].append([consensus_sequence, reference_name_consensus])
            else:
                consensus_dict[species] = [[consensus_sequence, reference_name_consensus]]

for species, data in consensus_dict.items():
    
    # if more than one entry, they will be aligned and a further consensus will be obtained later on
    # if not, output the sequence as the definitive consensus for the spp

    if len(data) > 1:

        destiny_directory_name = species.replace(" ","_") + "_consensus_directory"
        
        os.mkdir(destiny_directory_name, 0o777)
        
        for item in data:
            consensus_sequence = item[0]
            origin = os.path.realpath(consensus_sequence)
            destiny = f"{destiny_directory_name}/{consensus_sequence}"
            os.symlink(origin, destiny)

         
    else:
        consensus_sequence = data[0][0]

        origin = os.path.realpath(consensus_sequence)
        species_for_filename = species.replace(" ","_")
        destiny = f"./{species_for_filename}_consensus_sequence.fa"
        os.symlink(origin, destiny)





