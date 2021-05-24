#!/usr/bin/env python


import sys
import urllib.request
import argparse
import os

import random
random.seed(5)
############################################
####### 1. Parameters for the script #######
############################################

parser = argparse.ArgumentParser(description="Download assemblies from the NCBI refseq, already named for their immediate use in nf-core pikavirus. An internet network is required for this script to work.")

parser.add_argument("--only_ref_gen", type=bool, metavar="bool", default=False, dest="onlyref" , help="Download only those assemblies with the \"reference genome\" category (Default: False).")
parser.add_argument("--strains", type=bool , default=False,metavar="bool", dest="strains", help="Look for strains in the input taxids.")
parser.add_argument("-group", default="all", dest="group", choices=["virus","bacteria","fungi","all"], help="The group of assemblies from which download. Available: \"all\", \"virus\", \"bacteria\", \"fungi\" (Default: all).")
parser.add_argument("--only_sheet", type= bool,metavar="bool", default = False, dest="only_sheet", help="Download only the information sheet, in the proper order for pikavirus to work (Default: False).")
parser.add_argument("--only_complete_genome", type =bool, metavar="bool", default = False, dest="only_complete_genome", help="Download only those assemblies corresponding to complete genomes (Default: False).")
parser.add_argument("--single_assembly_per_species_taxid", type= bool, metavar="bool", default = False, dest="single_assembly_per_spp_taxid", help="Dowload only one assembly for each species taxid, reference-genome and complete genomes if able (Default: False).")
parser.add_argument("--single_assembly_per_subspecies_taxid", type= bool, metavar="bool", default = False, dest="single_assembly_per_strain_taxid", help="Dowload only one assembly for each subspecies/strain taxid, reference-genome and complete genomes if able (Default: False).")
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

    assembly_url = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/{group}/assembly_summary.txt"

    with urllib.request.urlopen(assembly_url) as response:
        assembly_file = response.read().decode("utf-8").split("\n")

    assembly_data = [line.split("\t") for line in assembly_file if not line.startswith("#")]
  
    filtered_assembly_data = [line for line in assembly_data if len(line)==22]

    # if activated, choose only complete genomes references
    if args.only_complete_genome == True:
        filtered_assembly_data = [line for line in filtered_assembly_data if line[11]=="Complete Genome"]

    # if activated, choose only reference genomes
    if args.onlyref == True:
        filtered_assembly_data = [line for line in filtered_assembly_data if line[4]=="reference genome"]

    # if activated, choose only one assembly per taxid (reference genomes and complete genomes are priority)
    if args.single_assembly_per_spp_taxid == True or args.single_assembly_per_strain_taxid == True:
        
        assembly_per_species_taxid = {}
        
        # one dict entry for each taxid, containing all assemblies refering to it
        # keys can be spp taxid or strain/subspp taxid

        for line in filtered_assembly_data:
            if args.single_assembly_per_spp_taxid == True:
                if line[5] not in assembly_per_species_taxid.keys():
                    assembly_per_species_taxid[line[5]] = [line]
                else:
                    assembly_per_species_taxid[line[5]].append(line)

            elif args.single_assembly_per_strain_taxid == True:
                if args.single_assembly_per_spp_taxid == True and args.single_assembly_per_strain_taxid == True:
                    print(f"If single_assembly_per_species_taxid and single_assembly_per_subspecies_taxid are both set to True, only subspecies/strains assemblies will be retrieved.")
                
                if line[6] not in assembly_per_species_taxid.keys():
                    assembly_per_species_taxid[line[6]] = [line]
                else:
                    assembly_per_species_taxid[line[6]].append(line)

        filtered_assembly_data = []

        # give a score to each assembly
        for _, data in assembly_per_species_taxid.items():
            
            chosen_list = []
            max_score = 0

            for datapiece in data:
                score = 0

                if datapiece[4] == "reference genome":
                    score += 6
                if datapiece[4] == "representative genome":
                    score += 4

                if datapiece[11] == "Complete Genome":
                    score += 5
                if datapiece[11] == "Chromosome":
                    score += 2

                if datapiece[13] == "Full":
                    score += 1

                if score > max_score:
                    max_score = score

                datapiece.append(score)

                chosen_list.append(datapiece)
            
            chosen_list = [item[0:-1] for item in chosen_list if item[-1] == max_score]
            chosen_assembly = random.choice(chosen_list)
            filtered_assembly_data.append(chosen_assembly)


    # remove repeated lines (avoid repeating downloads) and incomplete assemblies (should not happen at this point)
    filtered_assembly_data = [list(y) for y in set([tuple(x) for x in filtered_assembly_data])]
    assembly_quantity = len(filtered_assembly_data)
    
    # obtain the relevant information: 
    # assembly accession(0), subspecies taxid(1), species taxid(2), organism name(3), infraespecific name(4), assembly_level(5), ftp path(6)
    filtered_assembly_data = [[col[0],col[5],col[6],col[7],col[8],col[11],col[19]] for col in filtered_assembly_data]

    #############################################################################
    ####### 3. Generate the urls with an informative tsv for traceability #######
    #############################################################################

    for single_assembly in filtered_assembly_data:
        # generate the url for download
        ftp_path = single_assembly[-1]
        file_url = ftp_path.split("/")[-1]
        download_url = f"{ftp_path}/{file_url}_genomic.fna.gz"
        single_assembly.append(download_url)

    # generate the name for the file for download

        filename = f"{single_assembly[0]}.fna.gz"
        single_assembly.append(filename)

    # assembly accession(0)
    # subspecies taxid(1)
    # species taxid(2)
    # organism name(3)
    # infraespecific name(4)
    # assembly_level(5)
    # ftp path(6)
    # download_url(7)
    # filename(8)

    tsv_file = f"{group}_assemblies.tsv"

    with open(tsv_file,"w") as outfile:
        outfile.write(f"# Assemblies chosen for the group {group} \n")
        outfile.write(f"# Assembly_accession\tSpecies_taxid\tSubspecies_taxid\tScientific_name\tIntraespecific_name\tDownload_URL\tFile_name\tAssembly_level\n")
        for col in filtered_assembly_data:
            outfile.write(f"{col[0]}\t{col[2]}\t{col[1]}\t{col[3]}\t{col[4]}\t{col[7]}\t{col[8]}\t{col[5]}\n")

    print(f"{tsv_file} created successfully!")

    if args.only_sheet == False:

        ##############################################
        ####### 4. Download the assembly files #######
        ##############################################

        print(f"Starting download of {assembly_quantity} {group} assemblies!")
        destiny_folder = f"{group}_assemblies_for_pikavirus"
        
        if not os.path.exists(destiny_folder):
            os.mkdir(destiny_folder)

        for col in filtered_assembly_data:
            try:
                url = col[7]
                filename = col[8]
                location_filename = f"{destiny_folder}/{filename}"
                if os.path.exists(location_filename):
                    print(f"{filename} already found on destiny file.")
                else:
                    urllib.request.urlretrieve(url, location_filename)               

            except:
                print(f"Assembly {col[0]} from organism {col[3]} could not be retrieved.")

        print(f"Download of {group} assemblies complete!")

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