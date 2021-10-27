#!/usr/bin/env python

'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 1.0

CREATED: 30-7-2021

REVISED: 4-8-2021

DESCRIPTION: uses pandas and numpy to extract the coverage of the control genome
given to PikaVirus in the --control_sequence flag

INPUT (by order):
    -samplename: [str] name of the current sample (used for naming purposes)
    -coverage_file: [file/path] file containing the mapping of the control genome against the sample reads
    -idxstats: [file/path] file containing the idxstats of the control genome against the sample reads
    -flagstats: [file/path] file containing the flagstats of the control genome against the sample reads

OUTPUT:
    -{sample_name}_control_table.tsv: tsv file containing the coverage params of the control genome on the sample reads

USAGE:

    coverage_analysis_control.py samplename coverage_file idxstats flagstats

REQUIREMENTS: 
    -Pandas
    -Numpy
    -Python >= 3.6 

DISCLAIMER: this script has exclusively been developed for PikaVirus,
and therefore is not guaranteed to function properly in other settings. 
Despite this, feel free to use it at will.


TO DO: 

================================================================
END_OF_HEADER
================================================================
'''

# Imports
import sys
import os
import pandas as pd
import numpy as np

# args managent
sample_name = sys.argv[1]
coverage_file = sys.argv[2]
idxstats = sys.argv[3]
flagstats = sys.argv[4]

# Needed functions
def weighted_avg_and_std(df,values, weights):
    average = np.average(df[values], weights=df[weights])
    variance = np.average((df[values]-average)**2, weights=df[weights])
    
    return (average, variance**0.5)

def calculate_weighted_median(df, values, weights):
    cumsum = df[weights].cumsum()
    cutoff = df[weights].sum() * 0.5
    
    return df[cumsum >= cutoff][values].iloc[0]

# extract total number of reads
with open(flagstats) as total_reads:
    total_reads = total_reads.readlines()
    total_reads = [line.split("\t") for line in total_reads]

    print(total_reads)

    total_reads = int(total_reads[0][0])

# NC_001454.1	34214	12387	0
# *	0	0	1987613

with open(idxstats) as idxstats:
    idxstats = idxstats.readlines()
    idxstats = [line.split("\t") for line in idxstats]


# declare dict for final results
data = {"name":[],"covMean":[],"covSD":[],"covMin":[],"covMax":[],"covMedian":[],
        ">=x1":[],">=x10":[],">=x25":[],">=x50":[],">=x75":[],">=x100":[], "read_number":[], "total_reads":[]}

df = pd.read_csv(coverage_file,sep="\t",header=None)

df.columns=["gnm","covDepth","BasesAtThisCoverage","genomeLength","FracOnThisDepth"]

df["FracOnThisDepth_cumsum"] = df.groupby('gnm')['FracOnThisDepth'].transform(pd.Series.cumsum)
df["FracWithMoreDepth"] = 1 - df["FracOnThisDepth_cumsum"]
df["FracWithMoreDepth_percentage"] = df["FracWithMoreDepth"]*100

for name, df_grouped in df.groupby("gnm"):

    mean, covsd = weighted_avg_and_std(df_grouped,"covDepth","FracOnThisDepth")            
    minimum = min(df_grouped["covDepth"])
    maximum = max(df_grouped["covDepth"])
    median = calculate_weighted_median(df_grouped,"covDepth","FracOnThisDepth")

    # no full genome stats are shown, in case a multifasta is given as control sequence
    if name == "genome":
        continue
    else:
        gnm_name = f"Control sequence: {name}"

        for item in idxstats:
            if item[0] == name:
                sequence_reads = item[2]

    data["name"].append(gnm_name)
    data["covMean"].append(mean)
    data["covMin"].append(minimum)
    data["covMax"].append(maximum)
    data["covSD"].append(covsd)
    data["covMedian"].append(median)
    data[">=x1"].append(df_grouped.FracOnThisDepth[(df_grouped["covDepth"] >= 1)].sum())
    data[">=x10"].append(df_grouped.FracOnThisDepth[(df_grouped["covDepth"] >= 10)].sum())
    data[">=x25"].append(df_grouped.FracOnThisDepth[(df_grouped["covDepth"] >= 25)].sum())
    data[">=x50"].append(df_grouped.FracOnThisDepth[(df_grouped["covDepth"] >= 50)].sum())
    data[">=x75"].append(df_grouped.FracOnThisDepth[(df_grouped["covDepth"] >= 75)].sum())
    data[">=x100"].append(df_grouped.FracOnThisDepth[(df_grouped["covDepth"] >= 100)].sum())
    data["read_number"].append(sequence_reads)
    data["total_reads"].append(total_reads)


out = pd.DataFrame.from_dict(data)
out.to_csv(f"{sample_name}_control_table.tsv", sep="\t")


