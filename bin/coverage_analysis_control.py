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

REVISED: 30-7-2021

DESCRIPTION: 


INPUT (by order):

OUTPUT:

USAGE:

REQUIREMENTS:

DISCLAIMER: 


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


# Needed functions
def weighted_avg_and_std(df,values, weights):
    average = np.average(df[values], weights=df[weights])
    variance = np.average((df[values]-average)**2, weights=df[weights])
    
    return (average, variance**0.5)

def calculate_weighted_median(df, values, weights):
    cumsum = df[weights].cumsum()
    cutoff = df[weights].sum() * 0.5
    
    return df[cumsum >= cutoff][values].iloc[0]

# obtain name from filename



# declare dict for final results
data = {"name"=[],"covMean":[],"covSD":[],"covMin":[],"covMax":[],"covMedian":[],
        ">=x1":[],">=x10":[],">=x25":[],">=x50":[],">=x75":[],">=x100":[],"assembly":[]}

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

    file_prefix = coverage_file.replace("_coverage.txt").replace("_"," ")

    if name == "genome":
        gnm_name = f"control genome: {file_prefix}, full genome"
    else:
        gnm_name = f"control genome: {file_prefix}, sequence: {gnm}"

    data["gnm"].append(gnm_name)
    data["species"].append(species)
    data["subspecies"].append(subspecies)
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

out = pd.DataFrame.from_dict(data)
out.to_csv(f"{sample_name}_control_table.tsv", sep="\t")