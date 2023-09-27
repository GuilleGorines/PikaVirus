# /usr/bin/env python
'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII

AUTHOR: Guillermo J. Gorines Cordero

MAIL: gjgorines@isciii.es

VERSION: 0 

CREATED: 07-08-2023

DESCRIPTION: 
    Simple python-pandas script to transform csv files into excel files
    for delivery to normies

INPUT: 
    Text file, indicating the delimiter (e.g. if its a CSV, its a ","),
    and the desired name for the output

OUTPUT:

USAGE:
    python export_excel_from_csv.py 
        --input_file INPUT FILE
        --delimiter DELIMITER
        --output_filename NAME_OF_OUTPUT_FILE_NO_EXTENSION
        [--it_has_index]
        [--it_has_header]


REQUIREMENTS:
    -Python >= 3.6
    -Pandas


================================================================
END_OF_HEADER
================================================================
'''
import argparse
import pandas as pd
import sys


def parse_args(args=None):
    """
    Parse the args for the script
    """
    parser = argparse.ArgumentParser(
        description = "Transform csv or tsv file into an xlsx file, making sure that floats can be read as such"
    )

    parser.add_argument(
        "--input_file",
        dest = "input",
        required = True,
        type = str,
        help = "Input file, csv or tsv preferably"
    )

    parser.add_argument(
        "--delimiter",
        dest = "delimiter",
        default = ",",
        type = str,
        help = "Delimiter in the input file (default ',')"
    )

    parser.add_argument(
        "--output_filename",
        dest = "output",
        required = True,
        type = str,
        help = "Name desired for the output file (extension '.xlsx' will be added automatically)"
    )

    parser.add_argument(
        "--it_has_index",
        dest = "index",
        action = "store_true",
        type = bool,
        help = "whether or not the input file has index (if the rows are named)" 
    )

    parser.add_argument(
        "--it_has_header",
        dest = "header",
        action = "store_true",
        type = bool,
        help = "whether or not the input file has headers (if the columns are named)" 
    )

    return parser.parse_args(args)

def change_to_float(value):
    """
    HELPER FUNCTION for the convert_file function
    Change a value to a float if able, else leave it like that
    """
    try:
        result = float(value)
    except:
        result = value
    return result


def convert_file(
    infile,
    delimiter: str,
    header_in,
    index_in,
    outfile_prefix: str,
    ):
    """
    Load the file as a dataframe,
    Replace the data in the file as 
    """
    
    header = 0 if header_in else None
    index = 0 if index_in else None


    infile = pd.read_csv(
        infile,
        delimiter = delimiter,
        header = header,
        index_col = index
    )
    
    infile.apply(change_to_float, axis = 0)

    infile.to_excel(
        f"{outfile_prefix}.xlsx",
        header = header_in,
        index = index_in
    )

    return

def main(args=None):
    args = parse_args(args)
    convert_file(
        infile =  args.input,
        delimiter = args.delimiter,
        header_in = args.header,
        index_in = args.index,
        outfile_prefix = args.output,
        )

    return

if __name__ == "__main__":
    sys.exit(main())