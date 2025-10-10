''' EXTRACT read counts
Takes a folder of folders containing read count files from HTSEQ-count.
Prints read counts into one file.
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
parser.add_argument("infolder", type=str,
                    help="Input folder containing htseq-count output files (.txt)")
parser.add_argument("outfile", type=str,
                    help="Output file")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def list_files(infolder):
    '''Returns full paths of all .txt files in the folder sorted alphabetically'''
    files = [f for f in os.listdir(infolder) if f.endswith(".txt") and not f.endswith(".DS_Store")]
    files_sorted = sorted(files)  # sort alphabetically
    return [os.path.join(infolder, f) for f in files_sorted]
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    
    print("1. Extracting list of all genes expressed ...")
    files = list_files(args.infolder)
    print("Number of htseq count files found = " + str(len(files)))

    genenames = set()
    for file in files:
        with open(file, "r") as infile:
            for line in infile:
                line = line.rstrip()
                gene = line.split("\t")[0]
                # Adjust condition to match gene IDs in your files
                if not gene.startswith("__"):  
                    genenames.add(gene)
    print("Total number of genes with count data across all samples = " + str(len(genenames)))

    # Check that all files have same genes
    for file in files:
        sample_genes = []
        with open(file, "r") as infile:
            for line in infile:
                line = line.rstrip()
                gene = line.split("\t")[0]
                if not gene.startswith("__"):
                    sample_genes.append(gene)
        if len(sample_genes) != len(genenames):
            print("ERROR - gene set mismatch in file:", file)
            sys.exit(1)

    geneexpressiondict = defaultdict(list)
    header = []

    print("2. Extracting expression ...")
    for file in files:
        sample = os.path.basename(file).replace("_htseqcount.txt", "")
        header.append(sample)

        sampleexpressiondict = {}
        with open(file, "r") as infile:
            for line in infile:
                line = line.rstrip()
                gene = line.split("\t")[0]
                if not gene.startswith("__"):
                    counts = float(line.split("\t")[1])
                    sampleexpressiondict[gene] = counts

        print("No. of genes with count data in sample " + sample + ": " + str(len(sampleexpressiondict)))

        for gene in genenames:
            counts = sampleexpressiondict.get(gene, 0.0)  # default 0 if missing
            geneexpressiondict[gene].append(counts)

    print("Total number of genes for which read count info is extracted = " + str(len(geneexpressiondict)))

    print("3. Writing merged count output ...")
    with open(args.outfile, "w") as outfile:
        # Write header
        outfile.write("Geneid")
        for sample in header:
            outfile.write("\t" + sample)
        outfile.write("\n")

        # Write data
        count = 0
        for gene in sorted(geneexpressiondict.keys()):
            count += 1
            outfile.write(gene)
            for c in geneexpressiondict[gene]:
                outfile.write("\t" + str(c))
            outfile.write("\n")

    print("Number of genes written = " + str(count))


if __name__ == '__main__':
    main()