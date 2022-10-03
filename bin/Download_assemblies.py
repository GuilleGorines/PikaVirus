#!/usr/bin/env python
'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 1.1

CREATED: Exact date unknown (late 2020)

REVISED: 26-05-2021
         03-10-2022

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
        [-database]
        [--only_sheet]
        
        **SOON**
        [--only_ref_gen]
        [--only_complete_genome]
        [--single_assembly_per_species_taxid]
        [--single_assembly_per_subspecies_taxid]

REQUIREMENTS:
    -Python >= 3.6
    -Urllib

TO DO: 

================================================================
END_OF_HEADER
================================================================
'''
# General imports
import argparse
import os
import sys

# Specialiced imports
import urllib.request

# Variables for colors
info    = "\033[0;34m" + "INFO" + "\033[0m"
error   = "\033[0;31m" + "ERROR" + "\033[0m"
warning = "\033[1;33m" + "WARNING" + "\033[0m"
success = "\033[0;32m" + "SUCCESS" + "\033[0m"


# Local functions
def Get_db_assembly_list(database, group):
    """
    Get the assembly list 
    """
    reference = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/{database}/{group}/assembly_summary.txt"

    with urllib.request.urlopen(reference) as response:
        assembly_list = response.read().decode("utf-8").split("\n")

    assembly_data = [line.split("\t") for line in assembly_list if not line.startswith("#")]
    filtered_assembly_data = [line for line in assembly_data if len(line) == 23]

    """
    Anchor code, just in case the options "only complete genome" and "only ref" are accepted

    # if activated, choose only complete genomes references
    if args.only_complete_genome == True:
        filtered_assembly_data_refseq = [line for line in filtered_assembly_data_refseq if line[11]=="Complete Genome"]

    # if activated, choose only reference genomes
    if args.onlyref == True:
        filtered_assembly_data_refseq = [line for line in filtered_assembly_data_refseq if line[4]=="reference genome"]
    """
    filtered_assembly_data = [list(y) for y in set([tuple(x) for x in filtered_assembly_data])]
    assembly_quantity = len(filtered_assembly_data)
    
    # Obtain the relevant information: 
    # 0 - Assembly accession
    # 1 - Subspecies taxid
    # 2 - Species taxid
    # 3 - Organism name (scientific name)
    # 4 - Infraespecific name
    # 5 - Assembly_level
    # 6 - FTP path
    # 7 - repeated seqs the other db

    assemblies = [[col[0], col[5], col[6], col[7], col[8], col[11], col[19], col[17]] for col in filtered_assembly_data]

    return assemblies, assembly_quantity

def Join_refseq_gb(refseq_assembly, genbank_assembly):
    """
    Join the refseq and the genbank assembly files
    """
    redundant_assemblies = [col[7] for col in refseq_assembly]
    filtered_data_genbank = [line for line in genbank_assembly if line[0] not in redundant_assemblies]

    merged_assembly_data = refseq_assembly + filtered_data_genbank

    non_redundant_genbank_number = len(filtered_data_genbank)
    total_number = len(merged_assembly_data)

    return merged_assembly_data, total_number, non_redundant_genbank_number

def Write_assembly_data(assembly_data, group):
    """

    """
    with open(f"{}_assemblies.tsv","w") as outfile:
        outfile.write(f"# Assemblies chosen for the group {group} \n")
        outfile.write(f"# Assembly_accession\tSpecies_taxid\tSubspecies_taxid\tScientific_name\tIntraespecific_name\tDownload_URL\tFile_name\tAssembly_level\n")
        for col in assembly_data:
            outfile.write(f"{col[0]}\t{col[2]}\t{col[1]}\t{col[3]}\t{col[4]}\t{col[6]}\t{col[7]}\t{col[5]}\n")

    return

def Download_assemblies(assembly_data, group):
    """

    """
    # Check destiny folder
    destiny_folder = f"{group}_assemblies_for_pikavirus"
    if not os.path.exists(destiny_folder):
        os.mkdir(destiny_folder)

    for single_assembly in assembly_data:

        # Name for the download
        filename = f"{single_assembly[0]}.fna.gz" 
        location_filename = f"{destiny_folder}/{filename}"
        
        # URL name
        file_url = single_assembly[6].split("/")[-1]
        download_url = f"{ftp_path}/{file_url}_genomic.fna.gz"

        # Counter to see how many were correctly downloaded
        successful_downloads = 0
        unsuccessful_downloads = 0

        try:
            if os.path.exists(location_filename):
                print(f"{warning}: {filename} already found on destiny directory. Download skipped.")
            else:
                urllib.request.urlretrieve(url, location_filename)
                successful_downloads += 1     
        except:
            print(f"{error}: Assembly {single_assembly[0]} from organsim {single_assembly[3]} could not be retrieved. URL: {download_url}")
            unsuccessful_downloads += 1

        print(f"{success}: Download of {group} assemblies complete!")
        print(f"{info}: {successful_downloads} {group} assemblies downloaded successfully")
        print(f"{info}: {unsuccessful_downloads} {group} assemblies could not be downloaded")

    return 

############################################
####### 1. Parameters for the script #######
############################################

parser = argparse.ArgumentParser(
    description="Download assemblies from the NCBI refseq, already named for their immediate use in nf-core pikavirus. An internet network is required for this script to work."
    )

parser.add_argument(
    "-group",
    default="all", 
    dest="group", 
    choices=["virus","bacteria","fungi","all"],
    help="The group of assemblies from which download. Available: \"all\", \"virus\", \"bacteria\", \"fungi\" (Default: all)."
    )

parser.add_argument(
    "-database",
    default="all",
    dest="database",
    choices = ["refseq", "genbank", "all"],
    help="The database to download the assemblies from. Available: \"refseq\", \"genbank\",\"all\" (Default: all)."
    )

parser.add_argument(
    "--only_sheet",
    action='store_true',
    default = False,
    dest="only_sheet",
    help="Download only the information sheet, in the proper order for pikavirus to work (Default: False)."
    )
    
# parser.add_argument("--only_ref_gen", action='store_true', default=False, dest="onlyref" , help="Download only those assemblies with the \"reference genome\" category (Default: False).")
# parser.add_argument("--only_complete_genome", action='store_true', default = False, dest="only_complete_genome", help="Download only those assemblies corresponding to complete genomes (Default: False).")
# parser.add_argument("--single_assembly_per_species_taxid", action='store_true', default = False, dest="single_assembly_per_spp_taxid", help="Dowload only one assembly for each species taxid, reference-genome and complete genomes if able (Default: False).")
# parser.add_argument("--single_assembly_per_subspecies_taxid", action='store_true', default = False, dest="single_assembly_per_strain_taxid", help="Dowload only one assembly for each subspecies/strain taxid, reference-genome and complete genomes if able (Default: False).")
args = parser.parse_args()

#########################################################
####### 2. Get the refseq assembly data from NCBI #######
#########################################################

groups_to_download = ["virus","bacteria","fungi"] if args.group.lower() == "all" else [args.group.lower()]

# for each group

for group in groups_to_download:
    if group == "virus":
        group == "viral"
    
    if args.database == "all":
        assembly_list_refseq, refseq_number = Get_db_assembly_list("refseq", group)
        assembly_list_genbank, genbank_number = Get_db_assembly_list("genbank", group)
        merged_assembly_list, total_number, non_redundant_genbank_number  = Join_refseq_gb(assembly_list_refseq, assembly_list_genbank)
        Write_assembly_data(merged_assembly_list, group)

        print(f"{info}: {refseq_number} total RefSeq assemblies were found")
        print(f"{info}: {genbank_number} total GenBank assemblies were found")

        if args.only_sheet:
            print(f"{info}: {total_number} assemblies were included in the assembly sheet: ({refseq_number} from RefSeq; {non_redundant_genbank_number} from GenBank)")
            print(f"{success}: All done!")
            sys.exit(0)

        print(f"{info}: {total_number} assemblies will be downloaded ({refseq_number} RefSeq; {non_redundant_genbank_number} GenBank)")
        Download_assemblies(merged_assembly_list, group)

    else:
        assembly_list, number = Get_db_assembly_list(args.database, group)
        Write_assembly_data(assembly_list, group)

        nice_name = "RefSeq" if args.database == "refseq" else "GenBank"

        print(f"{info}: {number} total assemblies were found in {nice_name}")

        if args.only_sheet:
            print("All done!")
            sys.exit(0)

        print(f"{info}: {number} {nice_name} assemblies will be downloaded")
        Download_assemblies(merged_assembly_list, group)
    
    print(f"{success}: All done!")
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

# UNUSED CODE

# if activated, choose only one assembly per taxid (reference genomes and complete genomes are priority)

#if args.single_assembly_per_spp_taxid == True or args.single_assembly_per_strain_taxid == True:
        
#        assembly_per_species_taxid = {}
        
        # one dict entry for each taxid, containing all assemblies refering to it
        # keys can be spp taxid or strain/subspp taxid

#        for line in filtered_assembly_data:
#            if args.single_assembly_per_spp_taxid == True:
#                if line[5] not in assembly_per_species_taxid.keys():
#                    assembly_per_species_taxid[line[5]] = [line]
#                else:
#                    assembly_per_species_taxid[line[5]].append(line)

#            elif args.single_assembly_per_strain_taxid == True:
#                if args.single_assembly_per_spp_taxid == True and args.single_assembly_per_strain_taxid == True:
#                    print(f"If single_assembly_per_species_taxid and single_assembly_per_subspecies_taxid are both set to True, only subspecies/strains assemblies will be retrieved.")
                
#                if line[6] not in assembly_per_species_taxid.keys():
#                    assembly_per_species_taxid[line[6]] = [line]
#                else:
#                    assembly_per_species_taxid[line[6]].append(line)

#        filtered_assembly_data = []

        # give a score to each assembly
#        for _, data in assembly_per_species_taxid.items():
            
#            chosen_list = []
#            max_score = 0

#            for datapiece in data:
#                score = 0

#                if datapiece[4] == "reference genome":
#                    score += 6
#                if datapiece[4] == "representative genome":
#                    score += 4
#
#                if datapiece[11] == "Complete Genome":
#                    score += 5
#                if datapiece[11] == "Chromosome":
#                    score += 2

#                if datapiece[13] == "Full":
#                    score += 1

#                if score > max_score:
#                    max_score = score

#                datapiece.append(score)

#                chosen_list.append(datapiece)
            
#            chosen_list = [item[0:-1] for item in chosen_list if item[-1] == max_score]
#            chosen_assembly = random.choice(chosen_list)
#            filtered_assembly_data.append(chosen_assembly)
