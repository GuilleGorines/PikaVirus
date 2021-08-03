#!/usr/bin/env python
'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 

CREATED: 1-8-2021

REVISED: 1-8-2021

DESCRIPTION: 
    In nf-core pikavirus, generates the result HTML for each sample
    *insert description, third person*

INPUT (by order):
    1. *first_input*
        *note on first input*
    2. *second_input*

OUTPUT:


OPTIONS:


USAGE:
    *script_name*.py [options]

REQUIREMENTS:
    -Python >= 3.6

DISCLAIMER:

TO DO: 

================================================================
END_OF_HEADER
================================================================
'''

# Necessary imports
import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

msa_file = sys.argv[1]

fasta_name = msa_file.replace("_msa.fasta","")
header = f">{fasta_name}\n"

outfile_name = f"{fasta_name}_consensus_sequence_from_msa.fasta"

consensus_threshold = 0.5

alignment = AlignIO.read(msa_file, 'fasta')
summary_align = AlignInfo.SummaryInfo(alignment)
consensus_sequence = summary_align.gap_consensus(consensus_threshold)

with open(outfile_name,"w") as outfile:
    outfile.write(header)
    outfile.write(str(consensus_sequence))