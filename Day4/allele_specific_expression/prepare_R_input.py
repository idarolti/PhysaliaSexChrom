#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("autosomes", type=str,
					help="")
parser.add_argument("sexchromosomes", type=str,
					help="")
parser.add_argument("output", type=str,
					help="")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	count = 0
	with open(args.output, "w") as out:
		out.write("count,category,ratio\n")
		with open(args.autosomes, "r") as auto:
			next(auto)
			for line in auto:
				count += 1
				line = line.rstrip()
				fraction = line.split(",")[-1]
				out.write(str(count))
				out.write(",Autosomes,")
				out.write(str(fraction))
				out.write("\n")
		with open(args.sexchromosomes, "r") as sexchr:
			next(sexchr)
			for line in sexchr:
				count += 1
				line = line.rstrip()
				fraction = line.split(",")[-1]
				out.write(str(count))
				out.write(",Sex chromosomes,")
				out.write(str(fraction))
				out.write("\n")

if __name__ == '__main__':
	main()
