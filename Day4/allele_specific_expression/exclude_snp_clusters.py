import argparse
import sys
import vcfpy
import os
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, help="A VCF file")
parser.add_argument("-l", "--read_length", type=int, help="Length of the RNA-Seq reads")
parser.add_argument("-m", "--mismatches", type=int, help="Allowed mismatches in the alignment")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

def find_snp_cluster(snp_positions, read_length, mismatches):
    exclude = []
    for i in range(len(snp_positions)):
        cluster = [snp_positions[i]]
        for x in range(i+1, len(snp_positions)):
            if snp_positions[x] > snp_positions[i] + read_length:
                if len(cluster) > mismatches:
                    exclude += cluster
                break
            else:
                cluster.append(snp_positions[x])
    return set(exclude)

def get_snp_positions(vcf_reader):
    snp_positions = defaultdict(list)
    for record in vcf_reader:
        snp_positions[record.CHROM].append(record.POS)
    return snp_positions

def get_snp_cluster(snp_positions, read_length, mismatches):
    snp_cluster = defaultdict(list)
    for chromosome in snp_positions:
        snps = snp_positions[chromosome]
        snp_cluster[chromosome] = find_snp_cluster(snps, read_length, mismatches)
    return snp_cluster

def filter_clusters(infile, outfile, read_length, mismatches):
    vcf_reader = vcfpy.Reader.from_path(infile)
    snp_positions = get_snp_positions(vcf_reader)
    snp_cluster = get_snp_cluster(snp_positions, read_length, mismatches)

    vcf_reader = vcfpy.Reader.from_path(infile)
    vcf_writer = vcfpy.Writer.from_path(outfile, vcf_reader.header)

    for record in vcf_reader:
        if record.POS not in snp_cluster[record.CHROM]:
            vcf_writer.write_record(record)
    vcf_writer.close()

def main():
    infile = args.input
    head, tail = os.path.split(infile)
    outfile = tail.split(".")[0] + "_noclusters.vcf"
    outfile = os.path.join(head, outfile)
    print("Filtering file: %s" % infile)
    print(outfile)
    filter_clusters(infile, outfile, args.read_length, args.mismatches)

if __name__ == "__main__":
    main()