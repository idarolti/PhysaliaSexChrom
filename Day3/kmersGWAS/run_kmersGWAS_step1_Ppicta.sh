#!/bin/bash

# Script for running kmersGWAS from the SexFindR pipeline (https://sexfindr.readthedocs.io/en/latest/Fst.html)
# Running steps until create of kmers.table

conda activate sexchr

#####################################################################################################

# 1. Organise data structure

WORK_DIR="...Day3/kmersGWAS/Ppicta"
cd $WORK_DIR

# Make list of directories

for i in `cat dirlist.txt`
  do
  echo -n "Making directory $i"
  mkdir ./$i
  touch $i/input_files.txt
  echo -n "../${i}_R1.fastq.gz" > ${i}/input_files.txt
  echo -n -e "\n../${i}_R2.fastq.gz" >> ${i}/input_files.txt
  echo
  echo "Done."
done

for i in `cat dirlist.txt`
  do
  count=`find -name "${i}*" | wc -l`
  if [ $count -gt 3 ] || [ $count -lt 3 ]
  then
  echo
  echo -n "ERROR: Not all files for ${i} are present, copying did not run correctly!"
  exit 1
  fi
done

###########################################################################################################

# Execute step 1: kmc

echo "Running step 1: kmc"

for dir in $(cat dirlist.txt)
  do
  cd $dir
  echo
  echo -e "$dir"
  echo
  kmc -t16 -k31 -ci2 @input_files.txt output_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2
  kmc -t16 -k31 -ci0 -b @input_files.txt output_kmc_all ./ 1> kmc_all.1 2> kmc_all.2
  cd $WORK_DIR
done

for i in $(cat dirlist.txt)
  do
  if [ ! -f $i/output_kmc_all.kmc_pre ]
  then
  echo
  echo "ERROR: The output file $i/output_kmc_all.kmc_pre does not exist, kmc did not run correctly!"
  exit 1
  fi
done

#################################################################################################

# Execute step 2: add strand information

echo "Running step 2: Add strand information"

for dir in $(cat dirlist.txt)
  do
  cd $dir
  echo
  echo -e "$dir"
  echo
  /share/pool/CompGenomVert/Applications/kmersGWAS/bin/kmers_add_strand_information -c output_kmc_canon -n output_kmc_all -k 31 -o kmers_with_strand
  echo "Done."
  # echo "beginning cleanup"
  # rm *.kmc*
  # echo "done cleanup"
  cd $WORK_DIR
done

for dir in $(cat dirlist.txt) 
  do
  if [ ! -f $dir/kmers_with_strand ]
  then
  echo
  echo "ERROR: The output directory $i/kmers_with_strand does not exist, kmc did not run correctly!"
  exit 1
  fi
done

##################################################################################################

# Generate individual kmers list file

echo "Generating kmers list file for step 3"

START=$(date +%s)

cd $WORK_DIR

ls $WORK_DIR | tail -n +2 | awk '{printf "${WORK_DIR}/%s/kmers_with_strand\t%s\n", $NF,$NF}' > kmers_list_paths_1.txt
grep -v .gz kmers_list_paths_1.txt > kmers_list_paths_2.txt
grep -f dirlist.txt kmers_list_paths_2.txt > kmers_list_paths_clean.txt
sed 's-${WORK_DIR}/--g' kmers_list_paths_clean.txt > kmers_list_paths_final.txt
rm kmers_list_paths_clean.txt
rm kmers_list_paths_1.txt
rm kmers_list_paths_2.txt

# Step 3.1: Combine kmers

echo "Running step 3.1: Combining kmers lists into one file"

cd $WORK_DIR

bin/list_kmers_found_in_multiple_samples -l kmers_list_paths_final.txt -k 31 --mac 2 -p 0.2 -o kmers_to_use
# combine list of kmers
# kmers must appear in at least 2 individuals (--mac 2)
# kmer appears in either canonized or non-canoinized in at least 20% of the individuals it appears in (-p 0.2)

if [ ! -f kmers_to_use ] || [ ! -s kmers_to_use ]
  then
  echo
  echo "ERROR: The output file kmers_to_use does not exist or is empty, combining did not run correctly!"
  exit 1
fi

# Creating kmers table

echo "Running step 3.2: Creating table of kmers presence/absence"

cd $WORK_DIR

bin/build_kmers_table -l kmers_list_paths_final.txt -k 31 -a kmers_to_use -o kmers_table

if [ ! -f kmers_table.table ] || [ ! -s kmers_table.table ]
  then
  echo
  echo "ERROR: The output file kmers_table.table does not exist or is empty, table was not generated correctly!"
  exit 1
fi

echo "Completed step 3."
