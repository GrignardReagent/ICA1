# Make a directory called fastqc_out 
mkdir fastqc_out
# This script runs fastqc on each .gz file and saves the output in fastqc_out
fastqc -o fastqc_out/ -f fastq *.fq.gz
