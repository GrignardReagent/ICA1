##### STEP 3 PART 1 #####
echo -e "Welcome to STEP 3 PART 1, there is also a PART 2 to this."

##### PROCESS Step 1 #####
# Unzip all the .gz files (-d means decompress), we will need this later
echo -e "Unzipping all files..."
gzip -d *.gz
echo "All gz files decompressed and ready for processing"

##### PROCESS Step 2 #####
# Copy the Genome fasta file from /localdisk/data/BPSM/ICA1/Tcongo_genome/ to the current directory
# -r specifies copying the whole directory
cp -r /localdisk/data/BPSM/ICA1/Tcongo_genome/ .
 
cd ./Tcongo_genome
# Unzip the file
gzip -d TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz
# Move it to the previous directory
mv TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta ..
cd ..
# Clean up empty directory
rm -rf Tcongo_genome/

##### PROCESS Step 3 #####
# Run bowtie2 on TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta
# Index TriTrypDB-46_TcongolenseIL3000_2019 (No need to iterate this line…)
echo "Using _Genome.fasta as a reference sequence by running bowtie2-build"
# Building fasta_index, which we will use as a genome index for bowtie2
bowtie2-build TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta fasta_index
echo "Note: There are many index files generated, but we will use <fasta_index> as the genome index for the bowtie2 step!"
# In order for the for loop to work, we need to define two sets of files, the fqfile_1 for end 1 fq files
fqfile_1=*_1.fq
# Align reads to an indexed sequence using bowtie2
# Iterate through all fq files, feeding 2 variables into the bowtie2 function
for i in ${fqfile_1}
do
	# Extracting the sequence name by removing "_1.fq"
	fqname=${i%%_1.fq}
	# Construct the corresponding _2.fq filename within the for loop
	fqfile_2="${fqname}_2.fq"
	# Align reads to an indexed sequence using bowtie2
	# -1 and -2 were used as parameters as this is a pair-end sequencing, so we cannot use -U Tco-…fq; -x is the genome index fasta_index, -S is the output alignment in sam format
	echo "Alignment started... This can take a while."
	bowtie2 --no-unal -x fasta_index -1 ${i} -2 ${fqfile_2} -S "output_${fqname}.sam"
	echo "Alignment completed for sequence ${fqname}"
done

echo "End of STEP 3.1, now go to STEP 3.2 to use samtools."
