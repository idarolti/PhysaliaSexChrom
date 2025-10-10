#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("vcf", type=str)
parser.add_argument("positional_information", type=str)
parser.add_argument("out_autosomes", type=str)
parser.add_argument("out_sexchromosomes", type=str)
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".txt")]
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    scaffs_sexchromosomes = []
    scaffs_autosomes = []
    count_autosomes = 0
    count_sexchromosomes = 0

    with open(args.positional_information, "r") as pos:
        next(pos)
        for line in pos:
            scaffold = line.split(",")[0]
            chromo = line.split(",")[-3]
            if chromo == "chr8":
                scaffs_sexchromosomes.append(scaffold)
            else:
                scaffs_autosomes.append(scaffold)

    with open(args.out_autosomes, "w") as out_autosomes:
        with open(args.out_sexchromosomes, "w") as out_sexchromosomes:
            with open(args.vcf, "r") as vcf:
                for line in vcf:
                    if line.startswith("#"):
                        out_autosomes.write(line)
                        out_sexchromosomes.write(line)
                    else:
                        scaffold = line.split()[0]
                        if scaffold in scaffs_autosomes:
                            out_autosomes.write(line)
                            count_autosomes += 1
                        elif scaffold in scaffs_sexchromosomes:
                            out_sexchromosomes.write(line)
                            count_sexchromosomes += 1

    print("SNPs on autosomes = ", count_autosomes)
    print("SNPs on sex chromosomes = ", count_sexchromosomes)

if __name__ == '__main__':
    main()