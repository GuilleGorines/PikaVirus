#!/usr/bin/env python

'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 1.0

CREATED: Exact date unknown (early 2021)

REVISED: 26-5-2021

DESCRIPTION: 
    Extracts data from fastqc quality analysis

INPUT (by order):
    1. Sample prefix name
    2. Name of the directory containing the results
	3. "single_end" or "paired_end" (strings) to differentiate
	if single_end:
		4. fastqc result (txt) before trimming
		5. fastqc result (txt) after trimming
	if paired_end:
		4. fastqc result (txt) R1 before trimming
		5. fastqc result (txt) R2 before trimming
		6. fastqc result (txt) R1 after trimming
		7. fastqc result (txt) R2 after trimming

USAGE:
    extract_fastqc_data.py Sample_prefix Result_dir {"single_end","paired_end"} pre_result.txt post_result.txt

REQUIREMENTS:
    -Python >= 3.6

DISCLAIMER:

TO DO:
	-REVISAR EL PATH EN EL QUE DEBERIA SALIR, IGUAL NI DEBIERA INCLUIRSE E INCLUIRSE MÁS TARDE, CON LO QUE EL INPUT 2 SOBRARÁ

================================================================
END_OF_HEADER
================================================================
'''

import sys

# necessary function to extract data
def print_basic_report_data(report, post_pre,outfile):
	with open(report,"r") as infile:
		infile = infile.readlines()
		for line in infile:
			if line.startswith("Filename"):
				filename = line.replace("Filename\t","").replace("\n","")

			elif line.startswith("Total Sequences"):
				nseqs = line.replace("Total Sequences\t","").replace("\n","")

			elif line.startswith("Sequence length"):
				seqlen = line.replace("Sequence length\t","").replace("\n","")

			elif line.startswith("%GC"):
				gc_content = line.replace("%GC\t","").replace("\n","")

			html_file_name = pre_report.replace(".txt",".html")
			html_path =f"{result_dir}/raw_fastqc/{html_file_name}"

		outfile.write(f"{samplename},{single_end_statement},{post_pre},{filename},{seqlen},{nseqs},{gc_content},{html_path}\n")

		return

## Going sample by sample
## Sample name is supplied

samplename = sys.argv[1]
result_dir = sys.argv[2]
single_end = sys.argv[3]

if single_end == "True":
	pre_data = [sys.argv[4]]
	post_data= [sys.argv[5]]
	single_end_statement = "single_end"

else:
	pre_data = [sys.argv[4],sys.argv[5]]
	post_data= [sys.argv[6],sys.argv[7]]
	single_end_statement = "paired_end"


## Organize reports based on trimmed (post) or not yet (pre)

with open(f"{samplename}_quality.txt","w") as outfile:
	for pre_report in pre_data:
		print_basic_report_data(pre_report,"pre",outfile)

	for post_report in post_data:
		print_basic_report_data(post_report,"post",outfile)
