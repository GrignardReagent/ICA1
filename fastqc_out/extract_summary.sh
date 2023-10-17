#!/bin/bash
# Make sure we're in the correct directory
cd fastqc_out/
# Make a directory to receive our extracted files
mkdir extracted_summaries

# for loop iterating all zip files
for zip_file in *.zip
do

	# Use basename to get the dir name from the zip file name
	dir_name="$(basename ${zip_file} .zip)"

	# Define the file to extract
	summary="${dir_name}/summary.txt"

	# Extract the specified file from the zip file, 
	# -j is set so that unzip does not make new directories, but instead extract the summary.txt file to a directory called "extracted_summaries"
	# -d sets the destination directory for file extraction ("extracted_summaries"
	unzip -j ${zip_file} ${summary} -d extracted_summaries

	# Rename the extracted file based on the dir name 
	mv "extracted_summaries/summary.txt" "extracted_summaries/${dir_name}_summary.txt"
done

echo -e "All summary.txt files are extracted to directory extracted_summaries/"
