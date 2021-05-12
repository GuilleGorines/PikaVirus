#!/usr/bin/env python

# USAGE:
#
#  graphs_coverage.py Samplename coveragefiles
# 
# Calculates basic coverage statistics for coverage files provided. Samplename needed for file naming.
#
# This script has been developed exclusively for nf-core/pikavirus, and we cannot
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

species_data = [line.split("\t") for line in species_data]
species_data = [[line[3], line[4], line[6].strip("\n")] for line in species_data]


# Remove the extension of the file (so it matches the filename)
extensions = [".fna.gz",".fna"]

for item in species_data:
    for extension in extensions:
        if item[2].endswith(extension):
            item[2]=item[2].replace(extension,"")

for item in coverage_files:

    match_name_coverage = os.path.basename(item).split("_vs_")[0]

    for name in species_data:
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
    data["species"].append("".join(set(df_grouped["Species"])))
    data["subspecies"].append("".join(set(df_grouped["Subspecies"])))
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

    plt.savefig(f"{name}.pdf")
    plt.close()

newcov = pd.DataFrame.from_dict(data)
newcov.to_csv(f"{outfile_name}_table.csv")
