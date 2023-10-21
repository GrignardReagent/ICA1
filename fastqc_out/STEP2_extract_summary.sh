#!/bin/bash
##### PROCESS step 1 #####
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

echo "All summary.txt files are extracted to directory extracted_summaries/"
echo "They are all renamed <sequence_num_fastqc_summary.txt>"

##### PROCESS step 2 ##### 
echo "Now we create 3 separate files that contain the number of PASS/FAIL/WARN and the name of the sequence"

cd extracted_summaries
# Define where to write out the PASS/FAIL/WARN results
pass_report="pass_report.txt"
fail_report="fail_report.txt"
warn_report="warn_report.txt"
# Use grep and and a pipe to get the count of PASS/FAIL/WARN and save them in 3 separate files
# grep -c finds the keyword and counts it
grep -c "PASS" *.txt | awk '{print $1, $3}' >> ${pass_report}
echo "Passed entries have been saved in ${pass_report}"

grep -c "FAIL" *.txt | awk '{print $1, $3}' >> ${fail_report}
echo "Failed entries have been saved in ${fail_report}"

grep -c "WARN" *.txt | awk '{print $1, $3}' >> ${warn_report}
echo "Warning entries have been saved in ${warn_report}"

# View the fail counts for the sequences
head fail_report.txt
echo "Most of the sequences have 1 fail, and it is not reasonable to filter out data containing FAILs - we'd have nothing left to analyse going forward!"
# View the warn counts for the sequences
head warn_report.txt

echo "Nothing in particular flags up. No sequence will be filtered out. Continue!"
# Go back by 2 directory levels
cd ..
cd ..

echo "STEP 2 done, go to STEP3"
