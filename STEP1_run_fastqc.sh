##### STEP 1 #####
echo -e "Welcome to STEP 1"
# Make a directory called fastqc_out 
mkdir fastqc_out
# This script runs fastqc on each .gz file and saves the output in fastqc_out
fastqc -o fastqc_out/ -f fastq *.fq.gz
# Move STEP2_extract_summary.sh to fastqc_out
mv ../STEP2_extract_summary.sh ./fastqc_out
# Change directory to where all the output files are kept
cd fastqc_out
echo -e "STEP2_extract_summary.sh has been moved to the current working directory."
echo -e "End of STEP 1, please execute STEP2_extract_summary.sh next"
