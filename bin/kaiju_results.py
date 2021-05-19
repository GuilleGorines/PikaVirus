#!/usr/bin/env python

import sys
import pandas as pd
import plotly.colors
import plotly.offline
import plotly.express as px
import plotly.graph_objects as go

outfile_name = sys.argv[1]
file = sys.argv[2]


# Extraer los datos de id, tamaño del contig y coverage
def process_node_data(input_string, classified = False):
    output_list = []
    
    for item in input_string:
        data = item[1].split("_")
        node_id, node_length, node_coverage = data[1], data[3], data[5]
        
        if classified == True:
            taxid = item[2]
            score = item[3]
            
            identifiers = item[4].split(",")
            identifiers = [n for n in identifiers if len(n) > 0]         
            
            accession_numbers = item[5].split(",")
            accession_numbers = [n for n in accession_numbers if len(n) > 0]
            
            organism_name = item[7]
            status = "Classified"
            
        else:
            taxid = str()
            score = str()
            identifiers = str()
            accession_numbers=str()
            organism_name=str()
            status = "Unclassified"
            
        output_list.append([node_id, node_length, node_coverage, taxid, score, identifiers, accession_numbers, organism_name, status])
            
    return output_list

def plot_coincidences(classified_list):
    
    plot_dict = {}
    
    for item in classified_list:

        if item[8] == "Unclassified":
            if "Unclassified" not in plot_dict.keys():
                plot_dict["Unclassified"] = 1
            else:
                plot_dict["Unclassified"] += 1

        elif item[8] == "Classified":
            if item[7] not in plot_dict.keys():
                plot_dict[item[7]] = 1
            else:
                plot_dict[item[7]] += 1

    df = pd.DataFrame(plot_dict.items(), columns = ["Organism","Number of contigs"])
    
    df["Percentage (from total)"] = (df["Number of contigs"]*100)/df["Number of contigs"].sum()

    fig = px.pie(df,
                 values = "Number of contigs",
                 names="Organism",
                 hover_data=["Percentage (from total)"],
                 color_discrete_sequence=plotly.colors.sequential.Turbo,
                 title=f"{outfile_name}, kaiju identification result")

    plotly.offline.plot({"data": fig},
                        auto_open=False,
                        filename = f"{outfile_name}_kaiju_pieplot.html")

with open(file) as infile:
    infile = infile.readlines()
infile = [item for item in infile]
infile = [item.strip("\n").split("\t") for item in infile]

unclassified_list = [item for item in infile if item[0]=="U"]
classified_list = [item for item in infile if item[0]=="C"]

classified_treated = process_node_data(classified_list, classified=True)
unclassified_treated = process_node_data(unclassified_list)

all_treated = classified_treated.extend(unclassified_treated)

plot_coincidences(classified_treated)

with open(f"{outfile_name}_kaiju_result_classified.txt","w") as outfile:
    outfile.write("Node_ID\tNode_length\tNode_coverage\tMatch_taxid\tMatch_score\tIdentifiers\tAccession_number\tOrganism\n")
    for line in classified_treated:
        for subitem in line:
            if type(subitem) == list:
                outline = ",".join(subitem)
            else:
                outline=subitem
            outfile.write(outline)

            if line.index(subitem) == len(line)-1:
                outfile.write("\n")
            else:
                outfile.write("\t")
            
with open(f"{outfile_name}_kaiju_result_unclassified.txt","w") as outfile:
    outfile.write("Node_ID\tNode_length\tNode_coverage\n")
    for item in unclassified_treated:
        line = "\t".join(item[0:3])
        outfile.write(f"{line}\n")

