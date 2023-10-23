#!/bin/bash
##### STEP 5 #####
echo -e "Welcome to STEP 5, let's calculate the average of the counts per genes for each group"

# Our strategy here is to use an awk script within a for loop.
# Creating an array of file names
input_genecounts=(*"_genecounts.cov")
# Use a for loop, call each item in input_genecounts array each time, as "file"
for file in "${input_genecounts[@]}"
do
	# Set the group name 
	group_name="${file%%_genecounts.cov}"
	echo -e "Calculating average for ${group_name}"
	awk '{ # Variable to store the sum of field values
	sum = 0;
	# Variable to count the number of fields
	count=0;
	# Loop through fields from the 6th field, run as far as i reaches NF, i++ increasing the value of i by 1 after each loop
	for (i=6; i<= NF; i++){
		# Values in the columns are added up and stored in sum 
		sum += $i;
		# Count increases by 1 each time the loop is ran, giving us the number of sequences in the group
		count++;
}
# ? works like a if/else statement: avg is sum/count if count is non-zero, otherwise, avg is 0. This prevents dividing by 0.
avg = count > 0 ? sum / count : 0;
# $4 is the gene name, $5 is the gene description
print $4, $5, avg;
}' ${file} > "${group_name}_average.txt"
	echo -e "The result is stored in ${group_name}_average.txt"
done

echo -e "This is the end of STEP 5, please proceed to STEP 6"
