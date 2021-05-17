#!/usr/bin/env python

# USAGE:
#
#  graphs_coverage.py Samplename coveragefiles
# 
# Calculates basic coverage statistics for coverage files provided. Samplename needed for file naming.
#
# DISCLAIMER: This script has been developed exclusively for nf-core/pikavirus, and we cannot
# assure its functioning in any other context. However, feel free to use any part
# of it if desired.

# Imports
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Needed functions
def weighted_avg_and_std(df,values, weights):
    average = np.average(df[values], weights=df[weights])
    variance = np.average((df[values]-average)**2, weights=df[weights])
    
    return (average, variance**0.5)

def calculate_weighted_median(df, values, weights):
    cumsum = df[weights].cumsum()
    cutoff = df[weights].sum() / 2.
    
    return df[cumsum >= cutoff][values].iloc[0]


# args managent
outfile_name=sys.argv[1]
species_data=sys.argv[2]
coverage_files=sys.argv[3:]

with open(species_data) as species_data:
    species_data = species_data.readlines()

# Get headers and data
headers = [line.strip("\n").split("\t") for line in species_data if line.startswith("#")]
species_data = [line.strip("\n").split("\t") for line in species_data if not line.startswith("#")]

# Identify required columns through headers
file_headers = ["filename","file_name","file-name","file"]
species_name_headers = ["scientific_name","organism_name","organism","species_name","species"]
subspecies_name_headers = ["intraespecific_name","subspecies_name","strain","subspecies"]

for single_header in header:
    for item in single_header:
        if item.lower() in file_headers:
            file_column_index = single_header.index(file)

        elif item.lower() in species_name_headers:  
            species_column_index = single_header.index(file)

        elif item.lower() in subspecies_name_headers:
            subspecies_column_index = single_header.index(file)

# Exit with error status if one of the required groups is not identified
if not file_column_index or not species_column_index or not subspecies_column_index:
    if not file_column_index:
        print(f"No headers indicating \"File name\" were found.")
    if not species_column_index:
        print(f"No headers indicating \"Species name\" were found.")
    if not subspecies_column_index:
        print(f"No headers indicating \"Subspecies name\" were found.")

    print(f"Please consult the reference sheet format, sorry for the inconvenience!")
    sys.exit(1)



species_data = [[line[species_column_index], line[subspecies_column_index], line[file_column_index]] for line in species_data]

# Remove the extension of the file (so it matches the filename)
extensions = [".gz",".fna"]

species_data_noext = []

for item in species_data:
    for extension in extensions:
        filename_noext=item[2].replace(extension,"")
    species_data_noext.append([item[0],item[1],filename_noext])

for item in coverage_files:

    match_name_coverage = item.replace(".sam","").split("_vs_")[0]

    for name in species_data_noext:
        if name[2] == match_name_coverage:
            species = name[0]
            subspecies = name[1]

            with open(item,"r") as infile:
                infiledata = [line.strip("\n") for line in infile.readlines()]
                infiledata = [line.split("\t") for line in infiledata]
                                
            for line in infiledata:
                if line[0] == "genome":
                    line[0] = f"{species}_{subspecies}_genome"

            with open(item,"w") as outfile:
                for line in infiledata:
                    filedata ="\t".join(line)
                    outfile.write(f"{filedata}\t{species}\t{subspecies}\n")

dataframe_list = []

for filename in coverage_files:
    tmp_dataframe = pd.read_csv(filename,sep="\t",header=None)
    dataframe_list.append(tmp_dataframe)

df = pd.concat(dataframe_list)

df.columns=["gnm","covThreshold","fractionAtThisCoverage","genomeLength","diffFracBelowThreshold","Species","Subspecies"]

df["diffFracBelowThreshold_cumsum"] = df.groupby('gnm')['diffFracBelowThreshold'].transform(pd.Series.cumsum)
df["diffFracAboveThreshold"] = 1 - df["diffFracBelowThreshold_cumsum"]
df["diffFracAboveThreshold_percentage"] = df["diffFracAboveThreshold"]*100

data = {"gnm":[],"species":[],"subspecies":[],"covMean":[],"covMin":[],"covMax":[],"covSD":[],"covMedian":[],
        "x1-x4":[],"x5-x10":[],"x10-x19":[],">x20":[],"total":[]}

for name, df_grouped in df.groupby("gnm"):

    mean, covsd = weighted_avg_and_std(df_grouped,"covThreshold","diffFracBelowThreshold")
    
    if mean == 0:
        continue
    
    minimum = min(df_grouped["covThreshold"])
    maximum = max(df_grouped["covThreshold"])
    median = calculate_weighted_median(df_grouped,"covThreshold","diffFracBelowThreshold")

    data["gnm"].append(name)
    data["species"].append(species)
    data["subspecies"].append(subspecies)
    data["covMean"].append(mean)
    data["covMin"].append(minimum)
    data["covMax"].append(maximum)
    data["covSD"].append(covsd)
    data["covMedian"].append(median)
    
    y0=df_grouped.diffFracBelowThreshold[(df_grouped["covThreshold"] >= 1) & (df_grouped["covThreshold"] < 5)].sum()
    y1=df_grouped.diffFracBelowThreshold[(df_grouped["covThreshold"] >= 5) & (df_grouped["covThreshold"] < 10)].sum()
    y2=df_grouped.diffFracBelowThreshold[(df_grouped["covThreshold"] >= 10) & (df_grouped["covThreshold"] < 20)].sum()
    y3=df_grouped.diffFracBelowThreshold[(df_grouped["covThreshold"] >= 20)].sum()
    y4=y0+y1+y2+y3
    
    data["x1-x4"].append(y0)
    data["x5-x10"].append(y1)
    data["x10-x19"].append(y2)
    data[">x20"].append(y3)
    data["total"].append(y4)
    

    plt.figure() 
    df_grouped.plot.line(x="covThreshold",
                        y="diffFracAboveThreshold_percentage",
                        legend=None)
    plt.title(name)
    plt.xlabel("Coverage Threshold")
    plt.ylabel("% of reads above threshold")

    plt.savefig(f"{species}_{subspecies}_{name}.pdf")
    plt.close()

newcov = pd.DataFrame.from_dict(data)
newcov.to_csv(f"{outfile_name}_table.csv")
