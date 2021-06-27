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
import sys
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.offline

# Input management

species_data = sys.argv[1]
bedgraph_files = sys.argv[2:]

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

    # Take reference file name
    ref_name = bedgraph_file.split("_vs_")[0]

    # Remove extension
    for extension in extensions:
        ref_name = ref_name.replace(extension,"")

    for name in species_data_noext:
        if ref_name == name[2]:
            species = name[0]
            subspecies = name[1]

    if subspecies:
        fig_title = f"{species} {subspecies}"
    else:
        fig_title = f"{species}"

    # declare dict for data (inside: dicts for each subsequence)
    graph_dict = {}

    # open file and split lines
    with open(bedgraph_file) as infile:
        bedgraph = infile.readlines()
    bedgraph = [line.split() for line in bedgraph]

    # for line, add to dict 
    for name,begin,end,coverage in bedgraph:

        if name not in graph_dict.keys():
            # initiate dict for sequence
            graph_dict[name]={}

        for position in range(int(begin),int(end)):
            graph_dict[name][position] = int(coverage)

    # generate scaffold for the data
    fig = make_subplots(rows=len(graph_dict), 
                        cols=1,
                        x_title="Position",
                        y_title="Coverage depth")

    fig.update_layout(title_text = f"{fig_title} ({ref_name}), coverage depth by genome length")


    # position for the subplot
    position = 1

    for key, subdict in graph_dict.items():

        subfig = go.Scatter( x = list(subdict.keys()),
                             y = list(subdict.values()),
                             fill='tozeroy',
                             name = key)

        fig.append_trace(subfig,
                        row = position,
                        col = 1)
        
        fig_single = go.Figure()
        fig_single.add_trace(subfig)
        fig_single.update_layout(title=f"{fig_title}, sequence id: {key}, coverage depth by genome length")

        position += 1

        plotly.offline.plot({"data": fig_single},
                            auto_open = False,
                            filename = f"{key}_coverage_depth_by_pos.html")     


    plotly.offline.plot({"data": fig},
                    auto_open = False,
                    filename = f"{ref_name}_coverage_depth_by_pos.html")                  



    


