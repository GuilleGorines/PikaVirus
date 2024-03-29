#!/usr/bin/env python

'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 0.1

CREATED: 24-6-2021

REVISED: 

DESCRIPTION: 
    *insert description, third person*

INPUT (by order):
    1. *first_input*
        *note on first input*
    2. *second_input*

OUTPUT:


OPTIONS:


USAGE:
    generate_len_coverage_graph.py [options]

REQUIREMENTS:
    -Python >= 3.6

DISCLAIMER:

TO DO: 

================================================================
END_OF_HEADER
================================================================
'''

# Imports
import os
import sys
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.offline

# Input management

sample_name = sys.argv[1]
type_of_organism=sys.argv[2]
species_data = sys.argv[3]
bedgraph_files = sys.argv[4:]


# create directory to hold non-zero coverage files
destiny_folder = f"{sample_name}_valid_bedgraph_files_{type_of_organism}"
os.mkdir(destiny_folder, 0o777)

# open tsv with species data
with open(species_data) as species_data:
    species_data = species_data.readlines()

# Get headers and data
headers = [line.strip("\n").split("\t") for line in species_data if line.startswith("#")]
species_data = [line.strip("\n").split("\t") for line in species_data if not line.startswith("#")]

# Identify required columns through headers
file_headers = ["filename","file_name","file-name","file"]
species_name_headers = ["scientific_name","organism_name","organism","species_name","species"]
subspecies_name_headers = ["intraespecific_name","subspecies_name","strain","subspecies"]

for single_header in headers:
    for item in single_header:
        if item.lower() in file_headers:
            file_column_index = single_header.index(item)

        elif item.lower() in species_name_headers:  
            species_column_index = single_header.index(item)

        elif item.lower() in subspecies_name_headers:
            subspecies_column_index = single_header.index(item)

# Exit with error status if one of the required groups is not identified
if not file_column_index or not species_column_index or not subspecies_column_index:
    if not file_column_index:
        print(f"No headers indicating \"File name\" were found in the reference file.")
    if not species_column_index:
        print(f"No headers indicating \"Species name\" were found in the reference file.")
    if not subspecies_column_index:
        print(f"No headers indicating \"Subspecies name\" were found in the reference file.")

    print(f"Please consult the reference sheet format, sorry for the inconvenience!")
    sys.exit(1)

species_data = [[line[species_column_index], line[subspecies_column_index], line[file_column_index]] for line in species_data]

# Remove the extension of the file (so it matches the filename)
# species_data_noext contains species, subspecies, filename without extension

extensions = [".gz",".fna"]

species_data_noext = []

for item in species_data:
    filename_noext = item[2]
    for extension in extensions:
        filename_noext=filename_noext.replace(extension,"")

    species_data_noext.append([item[0],item[1],filename_noext])

for bedgraph_file in bedgraph_files:

    with open(bedgraph_file) as infile:
        bedgraph = infile.readlines()
        bedgraph = [line.split() for line in bedgraph]

    # discard if no coverage at all
    if len(bedgraph) == 1 and int(bedgraph[0][3]) == 0:
        continue

    # Take reference file name
    assembly_name = bedgraph_file.split("_vs_")[0]

    # Remove extension
    for extension in extensions:
        assembly_name = assembly_name.replace(extension,"")

    for name in species_data_noext:
        if assembly_name == name[2]:
            species = name[0]
            subspecies = name[1]

    if subspecies:
        spp = f"{species} {subspecies}"
    else:
        spp = f"{species}"

    # rename the origin file for posterior rescue
    origin = os.path.realpath(bedgraph_file)
    safe_spp = spp.replace(" ","_").replace("/","-")
    destiny = f"{destiny_folder}/{sample_name}_{safe_spp}_{assembly_name}_bedgraph.txt"
    os.symlink(origin, destiny)

    # declare dict for data (inside: dicts for each subsequence)
    graph_dict = {}

    # for line, add to dict 
    for name,begin,end,coverage in bedgraph:

        if name not in graph_dict.keys():
            # initiate dict for sequence
            graph_dict[name]={}

        for position in range(int(begin),int(end)):
            # covert to float before int so no value error raises
            graph_dict[name][position] = int(float(coverage))

    # generate scaffold for the data
    full_lenplot = make_subplots(rows=len(graph_dict), 
                                 cols=1,
                                 x_title="Position",
                                 y_title="Coverage depth")
    full_lenplot.update_layout(title_text = f"{sample_name}: {spp} ({assembly_name}), all sequences, coverage depth by genome length")


    full_lenplot_cap_500 = make_subplots(rows=len(graph_dict), 
                                         cols=1,
                                         x_title="Position",
                                         y_title="Coverage depth")
    full_lenplot_cap_500.update_layout(title_text = f"{sample_name}: {spp} ({assembly_name}), all sequences, coverage depth by genome length, capped to depth 500")

    

    figurename = f"{sample_name}: {spp} genome, depth distribution by single base"
    filename = f"{sample_name}_{safe_spp}_{assembly_name}_genome" 

    # position for the subplot
    position = 1

    for key, subdict in graph_dict.items():

        # generate the lenplot
        single_lenplot = go.Scatter( x = list(subdict.keys()),
                                     y = list(subdict.values()),
                                     fill='tozeroy',
                                     name = key)

        # add lenplot to the global data
        full_lenplot.append_trace(single_lenplot,
                                  row = position,
                                  col = 1)

        # to the cap 500 as well
        full_lenplot_cap_500.append_trace(single_lenplot,
                                          row = position,
                                          col = 1)
        full_lenplot_cap_500.update_yaxes(range=[0, 500], row=position, col=1)
                                  
        
        # add to a figure so it can be modified
        lenplot_single = go.Figure()
        lenplot_single.add_trace(single_lenplot)
        lenplot_single.update_layout(title=f"{sample_name}: {spp} ({assembly_name}) , sequence id: {key}, coverage depth by genome length")

        # generate html with the plot
        plotly.offline.plot({"data": lenplot_single},
                            auto_open = False,
                            filename = f"{sample_name}_{spp}_{assembly_name}_{key}_coverage_depth_by_pos.html".replace(" ","_").replace("/","-")) 

        # cap lenplot to 500 max in y axis
        lenplot_single_500 = go.Figure()
        lenplot_single_500.add_trace(single_lenplot)
        lenplot_single.update_layout(title=f"{sample_name}: {spp} ({assembly_name}) , sequence id: {key}, coverage depth by genome length")
        lenplot_single.update_yaxes(range=[0,500])

        plotly.offline.plot({"data": lenplot_single},
                    auto_open = False,
                    filename = f"{sample_name}_{spp}_{assembly_name}_{key}_coverage_depth_by_pos_capped_500.html".replace(" ","_").replace("/","-")) 




        position += 1
    

    plotly.offline.plot({"data": full_lenplot},
                         auto_open = False,
                         filename = f"{sample_name}_{spp}_{assembly_name}_full_coverage_depth_by_pos.html".replace(" ","_").replace("/","-"))

    plotly.offline.plot({"data": full_lenplot_cap_500},
                        auto_open = False,
                        filename = f"{sample_name}_{spp}_{assembly_name}_full_coverage_depth_by_pos_capped_500.html".replace(" ","_").replace("/","-"))                      



    


