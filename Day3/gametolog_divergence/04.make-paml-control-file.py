''' Prepare PAML
This script prepares sequences for PAML.
Takes each folder of gametolog sequences and removes all files apart from phylip file. 
'''
#==============================================================================
import argparse
import sys
import os
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="Infolder containing folders of gametologs")

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
            f = os.path.join(path, name)
            file_list.append(f)
    return file_list
#==============================================================================
#Main==========================================================================
#==============================================================================m
def main():
    count = 0 
    for orthogroup in list_folder(args.infolder):
        for f in list_files(orthogroup):
            if not f.endswith(".phy"):
                os.remove(f)
            else:
                count +=1                
                control_file = f[:-3]+"ctl"
                print(control_file)
                seq = f
                out_file = f[:-3]+"txt"

                try:
                    with open(control_file,"wb") as outfile:
                        outfile.write("      seqfile = "+seq+" * sequence data filename"+"\n")
                        outfile.write("      outfile = "+out_file+"           * main result file name"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen"+"\n")
                        outfile.write("      verbose = 1  * 0: concise; 1: detailed, 2: too much"+"\n")
                        outfile.write("        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below"+"\n")
                        outfile.write("        weighting = 0"+"\n")
                        outfile.write("        commonf3x4 = 0"+"\n")
                        outfile.write("        ndata = 1"+"\n")
                        outfile.write(""+"\n")

                except IOError:
                    print("File does not exit!")

    print("Number of folders modified =", count)

if __name__ == '__main__':
    main()
