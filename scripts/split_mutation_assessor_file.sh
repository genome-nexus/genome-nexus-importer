#!/bin/bash

# This file is used to split the mutationassessor4_for_genome_nexus.tsv file into 4 smaller files
# Because the original file is too large to be imported into MongoDB in one go
# And the script cannot extract .gz file, so we use gzip instead
# Only need to run this script once to generate the split files and upload to S3 bucket

# Unzip the original file
# Note: the original file is not included in the repo, it can be downloaded from the link: https://drive.google.com/file/d/1V6r65xJFF5fJ7b9JHwqkvCe8wWDrIBhd/view. 
# There is a copy stored in S3 bucket: ttps://genome-nexus-static-data.s3.amazonaws.com/mutationassessor4_for_genome_nexus.tsv.xz
echo "Unzipping mutationassessor4_for_genome_nexus.tsv.gz"
gunzip -c mutationassessor4_for_genome_nexus.tsv.gz > mutationassessor4_for_genome_nexus.tsv

# Split the file into 4 smaller files, excluding the header
echo "Splitting file into 4 parts"
HEADER=$(head -n 1 mutationassessor4_for_genome_nexus.tsv)
tail -n +2 mutationassessor4_for_genome_nexus.tsv | split -l $(($(wc -l < mutationassessor4_for_genome_nexus.tsv) / 4)) - mutationassessor_split_

# Add the header to each split file and gzip them
echo "Add header and gzip the split files"
for file in mutationassessor_split_*
do
  sed -i '' "1s/^/$HEADER\n/" $file
  mv "$file" "$file.tsv"
  gzip $file.tsv
done

# Clean up the original file
rm mutationassessor4_for_genome_nexus.tsv
echo "Splitting and compression completed."
