# This script copies the fastq files necessary for Step 1
##### STEP 0 #####
echo -e "Welcome to ICA1, this is ICA1 for BPSM. STEP 0 is required to copy in all files required for this exercise. Please do this if you want the rest of the codes to work. :)"
# Copying all data necessary for this ICA into the current working directory. -r is set for directories.
echo -e "Downloading the required data..."
cp -r /localdisk/data/BPSM/ICA1/fastq . 
echo -e "Moving STEP1_run_fastqc.sh to the fastq/ directory"
mv -f STEP1_run_fastqc.sh ./fastq
echo -e "Changing directory to fastq/"
cd fastq/ 
echo -e "STEP 0 complete, please execute STEP1_run_fastqc.sh next"
