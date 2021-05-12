#!/usr/bin/env python

import sys

# Usage:
# extract_reference_assemblies.py kraken2_report NCBI_assembly_summary.txt.
# The assembly summaries for different taxa are available in ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

krakenrep = sys.argv[1]
summary = sys.argv[2]
name_end = sys.argv[3]

with open(krakenrep) as krakenreport:
    krakenreport = [line.split("\t") for line in krakenreport.readlines()]
    idlist = []
    namelist = []
    for line in krakenreport:
        if line[3] == "S1":
            idlist.append(line[4])
            namelist.append(line[5].strip)
    for line in krakenreport:
        if line[3] == "S" and line[4] not in idlist:
            idlist.append(line[4])
            namelist.append(line[5].strip)
    
    if len(namelist) == 0:
        print(f"No species were identified in the kraken report.")
        sys.exit(2)    

# Report:
#   3: rank code (Unclass, Kingdom...)
#   4: taxID
#   5: scientific_name (indented)

with open(summary) as assembly_sum:
    assembly_sum = [line.split("\t") for line in assembly_sum.readlines() if not line.startswith("#")]

chosen_assemblies = []
for name in namelist:
    candidate = []
    try:
        candidate = [line for line in assembly_sum if line[8].split(":")[1] in name and line[7] == name]
    except:
        if len(candidate) == 0:
            candidate = [line for line in assembly_sum if line[7] == name and line[4] == "reference genome" or line[4] == "reference genome"]
            if len(candidate) == 0:
                candidate = [line for line in assembly_sum if line[7] == name]
        
        if len(candidate) > 0:
            chosen_assemblies.append(candidate[0])

assembly_sum = [[col[7],col[6],col[0],col[11],col[4],col[13],col[10],col[19]] for col in chosen_assemblies]
if len(assembly_sum) == 0:
    print(f"No {name_end} species were found in the kraken report.")
    sys.exit(2)

raw_assembly = [[col[0].replace(" ","_"),col[7]] for col in assembly_sum]

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

header=["Scientific_name", "Species_Taxonomic_ID", "Assembly_accession_chosen", "Assembly_level", "Refseq_category", "Representation", "Status", "Assembly_url"]
assembly_sum.insert(0,header)

namefile = f"chosen_assemblies_data_{name_end}.tsv"
with open(namefile, "w") as chosen_assemblies_file:
    for line in assembly_sum:
        for field in line:
            chosen_assemblies_file.write(str(field))
            if line.index(field) == len(line)-1:
                chosen_assemblies_file.write("\n")
            else:
                chosen_assemblies_file.write("\t")

url_file = f"url_download_{name_end}.sh"

with open(url_file,"w") as url_command_file:
    url_command_file.write(f'#/bin/bash \n')

    for sci_name,url in raw_assembly:
        filename = url.split("/")[-1]
        url = f"{url}/{filename}_genomic.fna.gz"
        url_command_file.write(f"wget --quiet -O {sci_name}_{name_end}.fna.gz {url}\n")
