#!/env/python


# 26-10-2021

# Imports
import sys

samplename = sys.argv[1]
list_of_sams = sys.argv[2:]

def digest_sam(list_of_sams):
    
    reads_dictionary = {}
    
    for sam in list_of_sams:
        dict_entry = sam.split(".fna.gz_vs_")[0]
        
        with open(sam, "r") as reads_list:
            reads_list = reads_list.readlines()
            # get first field, read name
            reads_list = [line.split("\t")[0] for line in reads_list]
            # remove duplicates
            reads_list = list(set(reads_list))
            
        reads_dictionary[dict_entry] = reads_list
        
    return reads_dictionary


def find_unique(reads_dictionary):   
    
    filename = f"{samplename}_reads_numbers.tsv"
    with open(filename,"w") as outfile:

        # header
        outfile.write("Identifier\tTotal reads mapping to this ref\tUnique reads\t% Unique\n")

        for key, readlist in reads_dictionary.items():

            # list with all reads
            all_else = list(reads_dictionary.values())

            # remove the first sample reads
            index_of_readlist = all_else.index(readlist)
            all_else.pop(index_of_readlist)

            # flatten list
            all_else=[item for bamfile in all_else for item in bamfile]   

            identifier = key.split("/")[-1]
            total_reads = len(readlist)
            
            # compare sets (readlist first)
            unique_reads = list(set(readlist).difference(set(all_else)))
            
            uniq_percentage = len(unique_reads)*100/total_reads

            outfile.write(f"{identifier}\t{len(readlist)}\t{len(unique_reads)}\t{uniq_percentage}\n")
            print(f"{key}, {total_reads} total reads, {len(unique_reads)} unique reads ({uniq_percentage} %).")
        
    return

reads_dictionary = digest_sam(list_of_sams)
find_unique(reads_dictionary)