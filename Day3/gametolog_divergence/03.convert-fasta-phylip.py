''' Convert FASTA FILE TO PHYLIP FORMAT
This script converts specified file in a folder of folders into phylip format.
All other files are removed. This step is necessary for Paml.
'''
#==============================================================================
import argparse
import sys
import os
from subprocess import Popen, list2cmdline
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="Infolder containing folders of gametologs")
parser.add_argument("uniq", type=str,
                    help="Unique identifier which file to convert ends with")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def exec_in_row(cmds):
    if not cmds:
        return  # empty list

    def done(p):
        return p.poll() is not None

    def success(p):
        return p.returncode == 0

    def fail():
        sys.exit(1)

    for task in cmds:
        p = Popen(task)
        p.wait()

    if done(p):
            if success(p):
                print("done!")
            else:
                fail()

def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if not f.endswith(".DS_Store")]

def list_files(current_dir):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if name.endswith(args.uniq):
                f = os.path.join(path, name)
                file_list.append(f)
            else:
                f = os.path.join(path, name)
                os.remove(f)
    return file_list
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    makeprankrun = []
    infiles = list_folder(args.infolder)
    for infile in infiles:
        for inf in list_files(infile):
            outpath = os.path.dirname(inf)
            outfile = outpath+"/"+os.path.dirname(inf).split("/")[-1].split("_")[1]
            makeprankrun.append(["prank", "-convert", "-d="+inf, 
                                               "-f=paml",
                                               "-o="+outfile,
                                               "-keep"
                                               ])

    # Passes the makeprankrun commands (a list of lists) to the exec_commands
    # function and executes them in parallel.
    exec_in_row(makeprankrun)
    
if __name__ == "__main__":
    main()
