#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import vcfpy
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
					help="A vcf file")
parser.add_argument("outfile", type=str,
					help="A vcf file filtered for triallelic SNPs, Ns and no genotype info")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	count_valid_SNPs = 0
	count_triallelic = 0
	count_N_SNPs = 0
	vcf_lines = 0
	bp = ["A","G","T","C"]

	vcf_reader = vcfpy.Reader.from_path(args.infile)
	vcf_writer = vcfpy.Writer.from_path(args.outfile, vcf_reader.header)

	for record in vcf_reader:
		vcf_lines += 1
		sample_DP = 0
		genotypes = []
		count_samples = 0
		# Remove INDEL info, only consider SNPS
		if record.is_snv(): 
			# Check INDEL filter has worked. REF should be 1bp
			if len(record.REF) == 1:
				# Check REF is not N
				if record.REF not in bp:
					count_N_SNPs += 1
				else:
					# Check not triallelic site
					if len(record.ALT) != 1:
						count_triallelic += 1
					else:
						for call in record.calls:  # changed attribute name
							count_samples += 1
							gt = call.data.get('GT')
							abq = call.data.get('ABQ', 0)
							rbq = call.data.get('RBQ', 0)

							# Ignore missing data (./.)
							if gt != "./.":
								# Ignore sites with quality < 20
								if gt == "0/1":
									if abq >= 20 and rbq >= 20:
										genotypes.append(gt)
								elif gt == "0/0":
									if rbq >= 20:
										genotypes.append(gt)
								elif gt == "1/1":
									if abq >= 20:
										genotypes.append(gt)
		if len(genotypes) == count_samples:
			count_valid_SNPs += 1
			vcf_writer.write_record(record)

	vcf_writer.close()

	print("No. of lines in vcf file =", vcf_lines)
	print("No. of valid SNPs =",count_valid_SNPs)
	print("No. of triallelic+ SNPs removed=", count_triallelic)
	print("No. of SNPS where REF=N removed=", count_N_SNPs)

if __name__ == '__main__':
	main()
