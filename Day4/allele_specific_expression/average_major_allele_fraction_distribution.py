
#==============================================================================
import argparse
import sys
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
					help="A vcf file")
parser.add_argument("out_density", type=str,
					help="outfile with info about the maj_allele_fraction for R plotting")
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
	with open(args.out_density, "w") as out_density:
		out_density.write("count,average_major_allele_ratio\n")
		with open(args.infile, "r") as vcf:
			for line in vcf:
				if line.startswith("#"):
					pass
				if line.startswith("scaff"):
					count_heterozygote = 0
					sum_fraction = 0
					line = line.rstrip()
					scaffold = line.split()[0]
					position = line.split()[1]
					samples = line.split()[9:12]
					for sample in samples:
						genotype = sample.split(":")[0]
						if genotype != "./.":
							RD = float(sample.split(":")[4])
							AD = float(sample.split(":")[5])
							if RD > AD:
								maj_allele_fraction = float(RD/(float(RD)+float(AD)))
							else:
								maj_allele_fraction = float(AD/(float(RD)+float(AD)))
							if genotype == "0/1":
								count_heterozygote += 1
								sum_fraction += maj_allele_fraction
					if count_heterozygote >= 1:
						count += 1
						average_fraction = sum_fraction/count_heterozygote
						out_density.write(str(count))
						out_density.write(",")
						out_density.write(str(average_fraction))
						out_density.write("\n")

if __name__ == '__main__':
	main()
