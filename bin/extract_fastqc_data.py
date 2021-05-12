#!/usr/bin/env python

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
