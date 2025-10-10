''' EXTRACT gene lengths
This script takes a file of read counts across samples extracted from HTseq-count. Extracts gene lengths and 
prints into one file.
'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("read_counts", type=str,
                    help="A txt file of read counts across samples")
parser.add_argument("gtf", type=str,
                    help="GTF file")
parser.add_argument("outfile", type=str,
                    help="An outfile")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def merge_coordinates(coordinates):
    sorted_by_lower_bound = sorted(coordinates, key=lambda tup: tup[0])
    merged = []

    for coord in sorted_by_lower_bound:
        if not merged:
            merged.append(coord)
        else:
            lastmergedcoord = merged[-1]
            lastmergedend = lastmergedcoord[1]
            lastmergedstart = lastmergedcoord[0]
            # test for intersection between current coords and merged coords:
            # we know via sorting that merged coord start <= current coord start
            currentstart = coord[0]
            currentend = coord[1]
            if currentstart <= lastmergedend:
                upper_bound = max(lastmergedend, currentend)
                merged[-1] = (lastmergedstart, upper_bound)  # replace by merged interval
            else:
                merged.append(coord)
    return merged

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    #get coordinates of genes in GTF file
    gene_coordinates = defaultdict(list)
    count = 0
    with open(args.gtf,"r") as infile:
        for line in infile:
            count += 1
            line = line.rstrip()
            if not line.startswith("#"):
                #want to ignore transcript as gives whole start-end 
                if not line.split("\t")[2] == "transcript":
                    scaffold = line.split("\t")[0]
                    gene = line.split("gene_id")[1].split('"')[1]
                    coordinate1 = float(line.split("\t")[3])
                    coordinate2 = float(line.split("\t")[4])
                    coord = [coordinate1, coordinate2]
                    gene_coordinates[gene].append(coord)
    print("Number of genes in GTF file = ", len(gene_coordinates))

    #merge coordinates to get start and length
    length_dict = {}
    for gene in gene_coordinates:
        coordinates = gene_coordinates[gene]
        if len(coordinates) == 1:
            start = float(coordinates[0][0])
            end = float(coordinates[0][1])
            length = float(end - start + 1)
            length_dict[gene] = length
        else:
            length = 0
            merged_coord = merge_coordinates(coordinates)
            start = float(merged_coord[0][0])
            end = float(coordinates[-1][1])
            for exon in merged_coord:
                tempstart = float(exon[0])
                tempend = float(exon[1])
                templength = float(tempend - tempstart + 1)
                length += templength
            length_dict[gene] = length
    print("Positional information extracted for genes =", len(length_dict))

    #print length for genes in same order as found in read count file
    count = 0
    with open(args.outfile, "w") as outfile:
        with open(args.read_counts, "r") as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith("MSTRG"):
                    gene = line.split("\t")[0]
                    length = length_dict[gene]
                    count += 1
                    outfile.write(gene)
                    outfile.write("\t")
                    outfile.write(str(length))
                    outfile.write("\n")
    print("Printing positional informational ...", count)
            

if __name__ == '__main__':
    main()
