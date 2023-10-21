#!/bin/bash
##### STEP 3.2 #####
##### PROCESS Step 1 #####
echo "You're now in STEP 3.2"
# Save all sam files into 1 variable so we can feed into for loop
output_file=output_*.sam

# Converting output.sam into a BAM alignment using samtools view
# BAM: Binary Alignment Map, which stores the same data as a compressed binary file (faster to work with)
for i in ${output_file}
do
	filename=${i%%.sam}
	# Convert output.sam into a BAM alignment using samtools view
	# -o specifies the name of the output file
	samtools view -o ${filename}.bam ${i}
	echo "Completed sam to bam for ${output_file}"
done

##### PROCESS Step 2 #####
# To allow us to create an index, we need to sort the alignment
bamfile=output_*.bam

# Using a for loop to feed variables into samtools
for i in ${bamfile}
do
	filename=${i%%.bam}
	# -O specifies the output format, -o specifies the sorted output file name
	# Followed by the name of our (unsorted) bam file
	samtools sort -O bam -o sorted_${i} ${i}
	echo "Completed sorting the alignment for ${filename} using samtool sort" 
done

##### PROCESS Step 3 ##### 
# The sorted BAM alignment file can now be indexed using samtools index
sortedbam=sorted_output_*.bam

# Using a for loop to feed variables into samtools
for i in ${sortedbam}
do
        filename=${i%%.bam}
        # Indexing the sorted bam files allows us to use bedtool multicov 
        samtools index ${i}
        echo "Completed indexing for ${filename} using samtool index"
done

