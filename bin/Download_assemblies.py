#!/usr/bin/env python
'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 1.0

CREATED: Exact date unknown (late 2020)

REVISED: 26-5-2021

DESCRIPTION: 
    Downloads the refseq assemblies needed by nf-core pikavirus (.fna.gz format),
    as well as the reference sheet for it to work.

OPTIONS:
    -group: download the assemblies for this group ("all","virus","bacteria","fungi")
    --only_ref_gen: download only reference genomes
    --only_sheet: do not download and only make the ref sheet
    --only_complete_genome: download only complete genomes
    --single_assembly_per_species_taxid: download only one ref per spp taxid
    --single_assembly_per_subspecies_taxid: download only one ref per subspp taxid

USAGE:
    Download_all_assemblies.py 
        [-group]
        [--only_ref_gen]
        [--only_sheet]
        [--only_complete_genome]
        [--single_assembly_per_species_taxid]
        [--single_assembly_per_subspecies_taxid]

REQUIREMENTS:
    -Python >= 3.6
    -Urllib
    -Biopython
    -gzip

TO DO: 

================================================================
END_OF_HEADER
================================================================
'''
from Bio import SeqIO
import gzip

import urllib.request

import sys
import argparse
import os

############################################
####### 1. Parameters for the script #######
############################################

parser = argparse.ArgumentParser(description="Download assemblies from the NCBI refseq, already named for their immediate use in nf-core pikavirus. An internet network is required for this script to work.")

parser.add_argument("-group", default="all", dest="group", choices=["virus","bacteria","fungi","all"], help="The group of assemblies from which download. Available: \"all\", \"virus\", \"bacteria\", \"fungi\" (Default: all).")
parser.add_argument("-database", default="all", dest="database", choices = ["refseq", "genbank", "all"], help="The database to download the assemblies from. Available: \"refseq\", \"genbank\",\"all\" (Default: all).")
args = parser.parse_args()

#########################################################
####### 2. Get the refseq assembly data from NCBI #######
#########################################################

# change virus for viral (thats how it figures in refseq)

if args.group.lower() == "all":
    groups_to_download = ["viral","bacteria","fungi"]

if args.group.lower() == "fungi":
    groups_to_download = ["fungi"]

if args.group.lower() == "bacteria":
    groups_to_download = ["bacteria"]

if args.group.lower() == "virus":
    groups_to_download = ["viral"]

# for each group
for group in groups_to_download:

    # access the assembly file from ncbi
    if args.database == "refseq" or args.database == "all":

        # generate url
        assembly_url_refseq = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/{group}/assembly_summary.txt"

        with urllib.request.urlopen(assembly_url_refseq) as response:
            assembly_file_refseq = response.read().decode("utf-8").split("\n")

        assembly_data_refseq = [line.split("\t") for line in assembly_file_refseq if not line.startswith("#")]

        filtered_assembly_data_refseq = [line for line in assembly_data_refseq if len(line)==23]

        # remove repeated lines (avoid repeating downloads) and incomplete assemblies (should not happen at this point)
        filtered_assembly_data_refseq = [list(y) for y in set([tuple(x) for x in filtered_assembly_data_refseq])]
        
        assembly_refseq_quantity = len(filtered_assembly_data_refseq)
        
        # if all databases, search for the entries that are equal in both (avoid redundance)
        if args.database == "all":
            genbank_in_refseq = [line[17] for line in filtered_assembly_data_refseq]

        # obtain the relevant information: 
        # assembly accession(0), species taxid(1), subspecies taxid(2), organism name(3), infraespecific name(4), ftp path(5)
        assemblies_refseq = [[col[0],col[5],col[6],col[7],col[8],col[19]] for col in filtered_assembly_data_refseq]

    # same process for genbank
    if args.database == "genebank" or args.database == "all":

        assembly_url_genbank = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/{group}/assembly_summary.txt"

        with urllib.request.urlopen(assembly_url_genbank) as response:
            assembly_file_genbank = response.read().decode("utf-8").split("\n")
        
        assembly_data_genbank = [line.split("\t") for line in assembly_file_genbank if not line.startswith("#")]
        filtered_assembly_data_genbank = [line for line in assembly_data_genbank if len(line) == 23]
    
        if args.database == "all":
            filtered_assembly_data_genbank = [line for line in filtered_assembly_data_genbank if line[0] not in genbank_in_refseq]
        
        filtered_assembly_data_genbank = [list(y) for y in set([tuple(x) for x in filtered_assembly_data_genbank])]

        assembly_genbank_quantity = len(filtered_assembly_data_genbank)
        
        # obtain the relevant information: 
        # assembly accession(0), species taxid(1), subspecies taxid(2), organism name(3), infraespecific name(4), ftp path(5)
        
        assemblies_genbank = [[col[0],col[6],col[5],col[7],col[8],col[19]] for col in filtered_assembly_data_genbank]
        
    if args.database == "all":
        final_assemblies = assemblies_refseq + assemblies_genbank
    elif args.database == "genbank":
        final_assemblies = assemblies_genbank
    elif args.database == "refseq":
        final_assemblies = assemblies_refseq

    # assembly accession (0)
    # species taxid (1)
    # subspecies taxid (2)
    # organism name (3)
    # infraespecific name (4)
    # ftp_path (5)

    final_assemblies_w_url = []
    
    for single_assembly in final_assemblies:
        # generate the url for download
        ftp_path = single_assembly[-1]
        file_url = ftp_path.split("/")[-1]
        download_url = f"{ftp_path}/{file_url}_genomic.fna.gz"

        single_assembly_url = single_assembly[:-1].append(download_url)
        final_assemblies_w_url.append(single_assembly_url)

    final_assemblies = final_assemblies_w_url

    del final_assemblies_w_url

    # assembly accession (0)
    # species taxid (1)
    # subspecies taxid (2)
    # organism name (3)
    # infraespecific name (4)
    # download_url (5)

    ##############################################
    ########  Download the assembly files ########
    ##############################################
    
    if args.database == "all":

        total_assemblies = assembly_refseq_quantity + assembly_genbank_quantity

        print(f"Starting download of {total_assemblies} {group} assemblies!\n \
                {assembly_refseq_quantity} assemblies from RefSeq\n \
                {assembly_genbank_quantity} assemblies from GenBank\n")

    elif args.database == "refseq":
            print(f"Starting download of {assembly_refseq_quantity} {group} assemblies from RefSeq!")

    elif args.database == "genbank":
            print(f"Starting download of {assembly_genbank_quantity} {group} assemblies from GenBank!")

    destiny_folder = f"{group}_assemblies_for_pikavirus"
    
    if not os.path.exists(destiny_folder):
        os.mkdir(destiny_folder)

    final_assemblies_w_subsequences = []

    for line in final_assemblies:
        try:
            url = line[5]

            filename = f"{line[0]}.fna.gz"

            location_filename = f"{destiny_folder}/{filename}"
            
            if os.path.exists(location_filename):
                print(f"{filename} already found on destiny file.")
            else:
                urllib.request.urlretrieve(url, location_filename)               

        except:
            print(f"Assembly {line[0]} from organism {line[3]} could not be retrieved. URL: {url}")

        try:
            assembly_subsequence_list = []
            assembly_length_list = []

            with gzip.open(location_filename,"rt") as assembly_gzip:
                for record in SeqIO.parse(assembly_gzip, "fasta"):
                    assembly_subsequence_list.append(record.id)
                    assembly_length_list.append(len(record.seq))
            
            assembly_w_subsequence = line.append(";".join(assembly_subsequence_list)).append(";".join(assembly_length_list))

            final_assemblies_w_subsequences.append(assembly_w_subsequence)

        except:
            print(f"Couldnt open {filename} to check subsequences.")
            pass
    
    print(f"Download of {group} assemblies complete!")      
    final_assemblies = final_assemblies_w_subsequences
    del final_assemblies_w_subsequences
    
    # assembly accession (0)
    # species taxid (1)
    # subspecies taxid (2)
    # organism name (3)
    # infraespecific name (4)
    # download_url (5)
    # Subsequences (6)
    # Subsequences length (7)

    with open(f"{group}_assemblies.tsv","w") as outfile:
        outfile.write(f"# Assemblies chosen for the group {group} \n")

        outfile.write(f"# Assembly_accession\tSpecies_taxid\tSubspecies_taxid\tScientific_name\tIntraespecific_name\tDownload_URL\tSubsequences_name\tSubsequences_length\n")
        
        for col in final_assemblies:
            outfile_line = "\t".join(col)
            outfile.write(f"{outfile_line}\n")

    print(f"{group} assembly data file created successfully!")
    
        

print("All done!")
sys.exit(0)

# Column 0: "assembly_accession"
# Column 1: "bioproject"
# Column 2: "biosample"
# Column 3: "wgs_master"
# Column 4: "refseq_category"
# Column 5: "taxid"
# Column 6: "(sub)species_taxid"
# Column 7: "organism_name"
# Column 8: "infraspecific_name"
# Column 9: "isolate"
# Column 10: "version_status"
# Column 11: "assembly_level"
# Column 12: "release_type"
# Column 13: "genome_rep"
# Column 14: "seq_rel_date"
# Column 15: "asm_name"
# Column 16: "submitter"
# Column 17: "gbrs_paired_asm"
# Column 18: "paired_asm_comp"
# Column 19: "ftp_path"
# Column 20: "excluded_from_refseq"
# Column 21: "relation_to_type_material"
# full info: ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt

