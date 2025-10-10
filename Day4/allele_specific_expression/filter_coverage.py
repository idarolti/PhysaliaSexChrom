#==============================================================================
import argparse
import sys
import vcfpy
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="A vcf file")
parser.add_argument("outfile", type=str,
                    help="A vcf file filtered for coverage")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Main==========================================================================
#==============================================================================

vcf_reader = vcfpy.Reader.from_path(args.infile)
vcf_writer = vcfpy.Writer.from_path(args.outfile, vcf_reader.header)

def genotype_passes(sample):
    gt = sample.data.get('GT')
    rd = sample.data.get('RD')
    ad = sample.data.get('AD')

    # If RD or AD missing, treat as pass
    if rd is None or ad is None:
        return True

    rd, ad = int(rd), int(ad)

    if rd >= ad:
        major_allele = rd
        minor_allele = ad
    else:
        major_allele = ad
        minor_allele = rd
    sum_alleles = rd + ad

    # Heterozygous check (0/1 or 1/0)
    if gt in ("0/1", "1/0"):
        return sum_alleles >= 15 and minor_allele >= 4
    # Homozygous check (0/0 or 1/1) -> pass if either RD or AD > 10
    elif gt == "0/0":
        return rd >= 15 and ad == 0
    elif gt == "1/1":
        return rd == 0 and ad >= 15
    else:
        return True

for record in vcf_reader:
    # Check all samples
    if all(genotype_passes(call) for call in record.calls):  # changed samples to calls
        vcf_writer.write_record(record)

vcf_writer.close()
