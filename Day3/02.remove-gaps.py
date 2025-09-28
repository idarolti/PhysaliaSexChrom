''' Removes GAPS IN ALIGNMENT
This script removes gaps in the alignment of each orthogroup. It also removes alignments
where the length after gap removal fails the cutoff threshold specified.
'''
#==============================================================================
import argparse
import sys
import os
import pickle
from Bio.Seq import Seq
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="An input folder containing folders of sequences of gametologs")
parser.add_argument("invalid_folder", type=str,
                    help="An output folder for invalid gametologs")
parser.add_argument("-cutoff", type=str,
                    help="bp")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if not f.endswith(".DS_Store")]

def list_files(current_dir):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if not name.endswith(".DS_Store"):
                if name.endswith("prank.best.fas"):
                    f = os.path.join(path, name)
                    file_list.append(f)
    return file_list

def read_fasta(source):
    SG = {}
    try:
        with open(source, "r") as file:
            for line in file.readlines():
                if line[0] == ">":
                    name = line[1:].rstrip().split()[0]
                    header = line[1:].rstrip()
                    SG[name] = ["", header]
                else:
                    SG[name][0] += line.rstrip()
            return SG
    except IOError:
        print("File does not exit!")

#==============================================================================
#Main==========================================================================
#==============================================================================m
def main():
    ammend = 0
    invalid_orthogroups = []
    for orthogroup in list_folder(args.infolder):
        for aln in list_files(orthogroup):
            ammend = ammend+1
            final_remove = {}
            fasta = read_fasta(aln)
            for s in fasta:
                for i, p in enumerate(fasta[s][0]):
                    if p == "-":
                        final_remove[i] = 0
            original_length = float(len(fasta[s][0]))
            final_length  = float(len(final_remove))
            if original_length - final_length < int(args.cutoff):
                invalid_orthogroups.append(orthogroup)
            else:
                outfilename = aln[:-8]+"gapsrm.fa"
                with open(outfilename,"w") as outfile:
                    for s in fasta:
                        seq = []
                        outfile.write(">"+s+"\n")
                        for i, fs in enumerate(fasta[s][0]):
                            if i not in final_remove:
                                seq.append(fs)
                                outfile.write(fs)
                        outfile.write("\n")

    print("%i orthogroups are invalid and were removed" % len(set(invalid_orthogroups)))
    print("Total number of orthogroups processed =", ammend)
    if not os.path.isdir(args.invalid_folder): os.mkdir(args.invalid_folder)
    for o in list(set(invalid_orthogroups)):
        import shutil
        shutil.move(o, args.invalid_folder)
    
if __name__ == '__main__':
    main()
