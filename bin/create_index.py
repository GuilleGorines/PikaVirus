#!/usr/bin/env python
'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: guillermo.gorines@urjc.es

VERSION: 

CREATED: 6-7-2021

REVISED: 

DESCRIPTION: 
    In nf-core pikavirus, generates the index for the run, with a link to the result of each sample

INPUT:


OUTPUT:
    -html file containing the 


OPTIONS:


USAGE:
    *script_name*.py [options]

REQUIREMENTS:
    -Python >= 3.6

DISCLAIMER:

TO DO: 

================================================================
END_OF_HEADER
================================================================
'''
import argparse
from datetime import date
from datetime import datetime

parser = argparse.ArgumentParser(description="Generates the result index HTML for nf-core pikavirus")

parser.add_argument("--quality-control", action='store_true', default=False, dest="quality_control", help="Has quality-control been activated?")
parser.add_argument("--control-removal", action='store_true', default=False, dest="control_removal", help="Has the control genome been removed?")

parser.add_argument("--kraken_scouting", dest="scouting", action='store_true', default=False, help="Has scouting been used?")

parser.add_argument("--virus", action='store_true', dest="virus",default=False, help="Was the coverage analysis performed for virus?")
parser.add_argument("--bacteria", action='store_true', dest="bacteria",default=False, help="Was the coverage analysis performed for bacteria?")
parser.add_argument("--fungi", action='store_true', dest="fungi",default=False, help="Was the coverage analysis performed for fungi?")

parser.add_argument("--translated-analysis",  action='store_true',default=False, dest="translated_analysis", help="Has translated analysis been used?")



parser.add_argument("--samplenames", nargs="+", required=True, dest="sample_list", help="Name of each sample, space separated")

args = parser.parse_args()

indexfile = "pikavirus_index.html"
date = date.today().strftime("%B %d, %Y")
hour = datetime.now().strftime("%H:%M")

if args.quality_control:
    quality_control="True"
else:
    quality_control="False"

if args.control_removal:
    control_removal="True"
else:
    control_removal="False"

if args.virus:
    virus="True"
    coverage_analysis="True"
else:
    virus="False"

if args.bacteria:
    bacteria="True"
    coverage_analysis="True"
else:
    bacteria="False"

if args.fungi:
    fungi="True"
    coverage_analysis="True"
else:
    fungi="False"

if args.translated_analysis:
    translated_analysis="True"
else:
    translated_analysis="False"

if args.scouting:
    scouting = "True"
else:
    scouting = "False"

with open(indexfile,"w") as outfile:

    # Start HTML
    outfile.write("<html>\n")

    # Head
    outfile.write("<head>\n \
                   <title>Pikavirus: result index</title>\n \
                   <meta charset=\"utf-8\">\n \
                   <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n \
                   <link href=\"https://cdn.jsdelivr.net/npm/bootstrap@5.0.1/dist/css/bootstrap.min.css\" rel=\"stylesheet\" integrity=\"sha384-+0n0xVW2eSR5OomGNYDnhzAbDsOXxcvSN1TPprVMTNDbiYZCxYbOOl7+AMvyTG2x\" crossorigin=\"anonymous\">\n \
                   <script src=\"https://cdn.jsdelivr.net/npm/bootstrap@5.0.1/dist/js/bootstrap.bundle.min.js\" integrity=\"sha384-gtEjrD/SeCtmISkJkNUaaKMoLD0//ElJ19smozuHV6z3Iehds+3Ulb9Bn9Plx0x4\" crossorigin=\"anonymous\"></script>\n \
                   <style>\n \
                   .icon { width: 32px; height: 32px;}\n \
                   td {border: 1px solid hsl(0, 0%, 87%); text-align: center; padding: 2%; width: 8%}\n \
                   .True { background-color:rgba(0,255,0,0.15); }\n \
                   .False { background-color:rgba(255, 0, 0, 0.15); }\n \
                   </style>\n \
                   </head>\n")
    
    # Start body
    outfile.write("<body>")

    # Title
    outfile.write("<div class=\"container-fluid\" style=\"background-color: #003B6F; color: white;\">\n \
                   <h1 style=\"padding: 60px; color:white; margin-top: 2%;\">Pikavirus: index</h1>\n \
                   </div>\n")


    wiki_link = ""
    
    # Navbar
    outfile.write(f"<nav class=\"navbar navbar-expand-md fixed-top navbar-dark bg-dark\">\n \
                   <div class=\"container-fluid\">\n \
                   <a class=\"navbar-brand\" style=\"color: white; font-size: 30px\">Result index </a>\n \
                   <button type=\"button\" class=\"navbar-toggler\" data-bs-toggle=\"collapse\" data-bs-target=\"#bs-example-navbar-collapse-1\" aria-expanded=\"False\">\n \
                   <span class=\"navbar-toggler-icon\">Toggle navigation</span>\n \
                   </button>\n \
                   <div class=\"collapse navbar-collapse\" id=\"navbar_buttons\">\n \
                   <ul class=\"navbar-nav me-auto mb-2 mb-md-0\">\n \
                   <li><a class=\"nav-link active\" href=\"{wiki_link}\" target=\"_blank\">How does this work?</a></li>\n \
                   </ul>\n \
                   <ul class=\"nav navbar-nav ml-auto\">\n \
                   <li><a class=\"nav-link active\" href=\"https://github.com/BU-ISCIII/PikaVirus\" target=\"_blank\"><img class=icon src=\"https://raw.githubusercontent.com/GuilleGorines/nf-core-pikavirus/a96b707ae03383d3cb727b935d348a0ed859f0c7/assets/github_logo.svg\"></a></li>\n \
                   <li><a class=\"nav-link active\" href=\"https://twitter.com/BUISCIII\" target=\"_blank\"><img class=icon src=\"https://raw.githubusercontent.com/GuilleGorines/nf-core-pikavirus/a96b707ae03383d3cb727b935d348a0ed859f0c7/assets/twitter_logo.svg\"></a></li>\n \
                   <li><a class=\"nav-link active\" href=\"https://www.facebook.com/pages/Escuela-Nacional-de-Sanidad-Isciii/203300096355772?fref=ts\" target=\"_blank\"><img class=icon src=\"https://raw.githubusercontent.com/GuilleGorines/nf-core-pikavirus/a96b707ae03383d3cb727b935d348a0ed859f0c7/assets/facebook_logo.svg\"></a></li>\n \
                   </ul>\n \
                   </div>\n \
                   </div>\n \
                   </nav>\n")

    # Divide page in two
    outfile.write(f"<div class=\"container-fluid\">\n \
                    <div class=\"row mb-2\">\n")

    global_multiqc_report = "multiqc_report.html"
    # First half: params used and time of creation
    outfile.write(f"<div class=\"container-fluid col-md-5\">\n \
                    <h2 style=\"padding: 3%\">Parameters</h2>\n \
                    <table id=\"params_table\">\n \
                    <tr>\n \
                    <td><a>Quality control</a></td>\n \
                    <td class={quality_control}><a>{quality_control}</a></td>\n \
                    </tr>\n \
                    <tr>\n \
                    <td><a>Control genome removal</a></td>\n \
                    <td class={control_removal}><a>{control_removal}</a></td>\n \
                    </tr>\n \
                    <tr>\n \
                    <td><a>Initial scouting</a></td>\n \
                    <td class={scouting}><a>{scouting}</a></td>\n \
                    </tr>\n \
                    <tr>\n \
                    <td><a>Coverage analysis</a></td>\n \
                    <td class={coverage_analysis}><a>{coverage_analysis}</a></td>\n \
                    </tr>\n \
                    <tr>\n \
                    <td><a>Virus</a></td>\n \
                    <td class={virus}><a>{virus}</a></td>\n \
                    </tr>\n \
                    <tr>\n \
                    <td><a>Bacteria</a></td>\n \
                    <td class={bacteria}><a>{bacteria}</a></td>\n \
                    </tr>\n \
                    <tr>\n \
                    <td><a>Fungi</a></td>\n \
                    <td class={fungi}><a>{fungi}</a></td>\n \
                    </tr>\n \
                    <tr>\n \
                    <td><a>Translated analysis</a></td>\n \
                    <td class={translated_analysis}><a>{translated_analysis}</a></td>\n \
                    </tr>\n \
                    </table>\n \
                    <div class=\"container-fluid\" style=\"padding: 3%;\">\n \
                    <p>This report was generated on {date}, at {hour}</p>\n \
                    </div>\n \
                    <a href=\"{global_multiqc_report}\" target=\"_blank\" class=\"list-group-item list-group-item-action\" style=\"text-align: center;\">View MultiQC file for all samples</a>\n \
                    </div>\n")


    # Second half: sample list
    outfile.write("<div class=\"container-fluid col-md-5\">\n \
                   <h2 style=\"padding: 3%\">Sample list</h2>\n \
                   <div class=\"list-group\">\n")

    for sample in sorted(args.sample_list):

        ref_result = f"{sample}_results.html"

        outfile.write(f"<a style=\"text-align: center;\" href={ref_result} target=\"_blank\" class=\"list-group-item list-group-item-action\">{sample}</a>")

    outfile.write("</div>\n \
                   </div>\n")

    # End of second half
    outfile.write("</div>\n \
                   </div>\n")

    # End of body
    outfile.write("</body>\n")

    # End of HTML
    outfile.write("</html>\n")
