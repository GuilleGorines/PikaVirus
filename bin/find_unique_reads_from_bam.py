#!/usr/bin/env python
'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: gjgorines@isciii.es

VERSION: 1

CREATED: 18-08-2023

REVISED: 28-08-2023

DESCRIPTION: 
    For a set of bam files (each for a reference), detect the number of reads that have 
    mapped, the number of reads that have not, and the number of reads that have ONLY
    mapped against said reference

USAGE:
    find_unique_reads.py [Samplename]
note: BAM files must be present for this script to work properly, and these must be
      named as PikaVirus does (this is, "{reference}_vs_{})

REQUIREMENTS:
    -pysam
    -Python >= 3.6

DISCLAIMER: this script has been created in the context of PikaVirus. Feel free to
use it or modify it at will
'''
import pysam
import glob
import sys

# read_dict:  read_name : [list of references for which the read mapped]
# ref_dict:   reference_name : [[mapped reads against this reference],
#                               [reads that mapped only to this reference],
#                               [reads that did not map against this reference]]

read_dict = dict()
ref_dict  = dict()

samplename = sys.argv[1]

# BAM files by pattern using glob
for bam_file in glob.glob("*.bam"):
    reference_name = bam_file.split("_vs_")[0].replace(".fna.gz","")
    
    # list 1: reads that mapped (mapped)
    # list 2: reads that mapped ONLY to this very reference (unique)
    # list 3: reads that did NOT map (unmapped)
    ref_dict[reference_name] = [list(), list(), list()]

    # Open BAM file with pysam
    bamfile_handle = pysam.AlignmentFile(bam_file, "rb")
    infile = [str(item) for item in list(bamfile_handle.fetch(until_eof=True))]

    # Remove header
    alignment_lines = [line.split("\t") for line in infile if not line.startswith("@")]

    unmapped_flag_values = ["69", "73", "89", "101", "117", "121", "133", "137", "153", "165", "181", "185"]

    for alignment in alignment_lines:
        # if read name (alignment[0]) is not a key in the read_dictionary, add it
        if alignment[0] not in read_dict.keys():
            read_dict[alignment[0]] = []
        
        # if FLAG value is different to any of the values, it means it mapped against this reference
        # source: https://samtools.github.io/hts-specs/SAMv1.pdf, page 7
        # https://www.samformat.info/sam-format-flag

        if str(alignment[1]) not in unmapped_flag_values:
            read_dict[alignment[0]].append(reference_name)

        # Add unmapped reads to the list
        # this means that for that reference
        # the read with that identifier did NOT map
        else:
            ref_dict[reference_name][2].append(alignment[0])
    # Close the bam file
    bamfile_handle.close()

# Calculate total number of reads AFTER all bam files have been processed
total_number_reads = len(read_dict.keys())

# parse the read_dict
# if a read only mapped one single reference, add it to the ref_dict
# first and second lists (mapped and unique mapped)

for read, references in read_dict.items():
    if len(set(references)) == 1:
        ref_dict[references[0]][0].append(read)
        ref_dict[references[0]][1].append(read)
    else:
        for single_reference in references:
            ref_dict[single_reference][0].append(read)

with open(f"{samplename}_mapping_balance.tsv", "w") as outfile_stats: 

    # Header for stats file
    outfile_stats.write("Samplename\tReference Name\tNumber of mapped reads\tNumber of unique reads\tNumber of unmapped reads\tTotal number of reads in the sample\n")

    for reference, reads in ref_dict.items():
        mapped_reads_number = len(set(reads[0]))
        unique_reads_number = len(set(reads[1]))
        unmapped_reads_number  = len(set(reads[2]))
        outfile_stats.write(f"{samplename}\t{reference}\t{mapped_reads_number}\t{unique_reads_number}\t{unmapped_reads_number}\t{total_number_reads}\n")

        # If we ever wanted to keep the mapped reads
        # all_mapped_reads = "\n".join(reads[0]) if len(reads[0]) > 1 else str(reads[0])
        # with open(f"{samplename}_{reference}_mapped_reads.txt", "w") as outfile_mapped:
        #    outfile_mapped.write(all_mapped_reads)

        # reads with a single incidence
        # If we ever wanted to keep the unique reads
        # unique_reads = "\n".join(reads[1]) if len(reads[1]) > 1 else str(reads[0])
        # with open(f"{samplename}_{reference}_unique_reads.txt", "w") as outfile_unique:
        #    outfile_unique.write(unique_reads)

        # reads that did not map 
        # If we ever wanted to keep unmapped reads
        # if len(reads[2]) == 0:
        #     unmapped_reads = "NaN"
        # else:
        #     unmapped_reads = "\n".join(reads[2]) if len(reads[2]) > 1 else str(reads[0])
        # with open(f"{samplename}_{reference}_unmapped_reads.txt", "w") as outfile_unmapped:
        #     outfile_unmapped.write(unmapped_reads)

sys.exit(0)
