import sys
import urllib.request
import argparse
import os


############################################
####### 1. Parameters for the script #######
############################################

parser = argparse.ArgumentParser(description="Download assemblies from the NCBI refseq, already named for their immediate use in nf-core pikavirus. An internet network is required for this script to work.")

parser.add_argument("--only_ref_gen", type=bool, metavar="bool", default=False, dest="onlyref" , help="Download only those assemblies with the \"reference genome\" category (Default: False).")
parser.add_argument("--strains", type=bool , default=False,metavar="bool", dest="strains", help="Look for strains in the input taxids.")
parser.add_argument("-group", default="all", dest="group", choices=["virus","bacteria","fungi","all"], help="The group of assemblies from which download. Available: \"virus\", \"bacteria\", \"fungi\" (Default: all).")
parser.add_argument("--only_sheet", type= bool,metavar="bool", default = False, dest="only_sheet", help="Download only the information sheet, in the proper order for pikavirus to work.")
parser.add_argument("--only_complete_genome", type =bool, metavar="bool", default = False, dest="only_complete_genome", help="Download only those assemblies corresponding to complete genomes.")
args = parser.parse_args()

#########################################################
####### 2. Get the refseq assembly data from NCBI #######
#########################################################

# 2.1: change virus for viral (thats how it figures in refseq)

if args.group.lower() == "all":
    groups_to_download = ["viral","bacteria","fungi"]

if args.group.lower() == "fungi":
    groups_to_download = ["fungi"]

if args.group.lower() == "bacteria":
    groups_to_download = ["bacteria"]

if args.group.lower() == "virus":
    groups_to_download = ["viral"]

# 2.2: access the assembly file

for group in groups_to_download:
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

# 3.3: remove repeated lines (avoid repeating downloads) and incomplete assemblies
    filtered_assembly_data = [list(y) for y in set([tuple(x) for x in filtered_assembly_data])]
    assembly_quantity = len(filtered_assembly_data)
    
# 3.4: obtain the relevant information: 
# assembly accession(0), subspecies taxid(1), species taxid(2), organism name(3), infraespecific name(4), assembly_level(5), ftp path(6)
    filtered_assembly_data = [[col[0],col[5],col[6],col[7],col[8],col[11],col[19]] for col in filtered_assembly_data]

    #############################################################################
    ####### 4. Generate the urls with an informative tsv for traceability #######
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
    ####### 5. Download the assembly files #######
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
# Column 6: "species_taxid"
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