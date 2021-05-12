#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    "nf-core/pikavirus": ["v_pipeline.txt", r"(\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
    "FastQC": ["v_fastqc.txt", r"FastQC v(\S+)"],
    "Fastp": ["v_fastp.txt", r"fastp\s(\S+)"],
    "Kraken2": ["v_kraken2.txt", r"Kraken\sversion\s(\S+)"],
    "Kaiju": ["v_kaiju.txt", r"Kaiju (\S+)"],
    "Bowtie2": ["v_bowtie2.txt", r"bowtie2-align-s\sversion\s(\S+)"],
    "Mash": ["v_mash.txt",r"Mash\sversion\s(\S+)"],
    "Samtools" : ["v_samtools.txt",r"samtools\s(\S+)"],
    "BedTools": ["v_bedtools.txt",r"bedtools\sv(\S+)"],
    "Spades": ["v_spades.txt",r"SPAdes\sv(\S+)"],
    "Quast": ["v_quast.txt",r"QUAST\sv(\S+)"]
}

results = OrderedDict()
results["nf-core/pikavirus"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results["FastQC"] = '<span style="color:#999999;">N/A</span>'
results["Fastp"] = '<span style="color:#999999;">N/A</span>'
results["Kraken2"] = '<span style="color:#999999;">N/A</span>'
results["Kaiju"] = '<span style="color:#999999;">N/A</span>'
results["Bowtie2"] = '<span style="color:#999999;">N/A</span>'
results["Mash"] = '<span style="color:#999999;">N/A</span>'
results["Samtools"] = '<span style="color:#999999;">N/A</span>'
results["BedTools"] = '<span style="color:#999999;">N/A</span>'
results["Spades"] = '<span style="color:#999999;">N/A</span>'
results["Quast"] = '<span style="color:#999999;">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/pikavirus Software Versions'
section_href: 'https://github.com/nf-core/pikavirus'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
