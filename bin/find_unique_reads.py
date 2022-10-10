
import sys
import glob


read_dict = dict()
ref_dict  = dict()

samplename = sys.argv[1]

# SAM files by pattern
for sam_file in glob.glob("*.sam"):

    reference_name = sam_file.split("_vs_")[0].replace(".fna.gz","")
    
    # list 1: reads that mapped (mapped)
    # list 2: reads that mapped ONLY to this very reference (unique)
    # list 3: reads that did NOT map (unmapped)
    ref_dict[reference_name] = [list(), list(), list()]

    with open(sam_file, "r") as infile:
        # remove header
        alignment_lines = [line.split("\t") for line in infile if not line.startswith("@")]

    for alignment in alignment_lines:

        if alignment[0] not in read_dict.keys():
            read_dict[alignment[0]] = []
        
        if str(alignment[1]) != "4":
            read_dict[alignment[0]].append(reference_name)
        # add unmapped reads to the list
        else:
            ref_dict[reference_name][2].append(alignment[0]) 

total_number_reads = len(read_dict.keys())

for read, references in read_dict.items():

    if len(set(references[0])) == 0:
        unmapped_reads.append(read)
    
    elif len(set(references)) == 1:
        ref_dict[references[0]][0].append(read)
        ref_dict[references[0]][1].append(read)
    else:
        for single_reference in references:
            ref_dict[single_reference][0].append(read)

with open(f"{samplename}_mapping_balance.tsv", "w") as outfile_stats: 

    # Header for stats file
    outfile_stats.write("Samplename\tReference Name\tNumber of mapped reads\tNumber of unique reads\tNumber of unmapped reads\tTotal number of reads in the sample\n")

    for reference, reads in ref_dict.items():
                
        # reads that mapped 
        all_mapped_reads        = "\n".join(reads[0]) if len(reads[0]) > 1 else str(reads[0])
        all_mapped_reads_number = len(reads[0])
        with open(f"{samplename}_{reference}_mapped_reads.txt", "w") as outfile_mapped:
            outfile_mapped.write(all_mapped_reads)

        # reads with a single incidence
        unique_reads            = "\n".join(reads[1]) if len(reads[1]) > 1 else str(reads[0])
        unique_reads_number     = len(reads[1])
        with open(f"{samplename}_{reference}_unique_reads.txt", "w") as outfile_unique:
            outfile_unique.write(unique_reads)

        # reads that did not map 
        if len(reads[2]) == 0:
            unmapped_reads = "NaN"
        else:
            unmapped_reads = "\n".join(reads[2]) if len(reads[2]) > 1 else str(reads[0])
        unmapped_reads_number  = len(reads[2])

        with open(f"{samplename}_{reference}_unmapped_reads.txt", "w") as outfile_unmapped:
            outfile_unmapped.write(unmapped_reads)

        outfile_stats.write(f"{samplename}\t{reference}\t{all_mapped_reads_number}\t{unique_reads_number}\t{unmapped_reads_number}\t{total_number_reads}\n")


sys.exit(0)