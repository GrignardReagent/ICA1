#!/bin/bash
##### STEP 4 #####
echo "You are now in STEP 4"
##### PROCESS Step 1 #####
echo -e "Fetching the required file"
cp /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed .

# We need to group the sequences from this step, so that we can simply use bedtool multicol instead of having to bother retrieving and adding stuff from one file to aother!
# Input file with data
input_file="Tco2.fqfiles"

# Sort the input data by relevant columns (SampleType, Replicate, Time, Treatment)
# This step is necessary for grouping the data effectively.
sorted_file="sorted_${input_file}"
# Using sort to specify sorting keys based on columns 2, 4 and 5, sorting priority: column 2, 4 then 5.
# -t$'\t' specifies the tab character as a field delimiter
awk 'NR>1' "${input_file}" | sort -t$'\t' -k2,2 -k4,4 -k5,5 > "${sorted_file}"

# Set an empty string variable so that the if loop can work
current_group=""
group_file="groups_Tco2.txt"
final_group_file="final_groups_Tco2.txt"
# -a assigns the words read to sequential indices of the array variable ARRAY, starting at 0, 
# -r makes sure that no backslashes are allowed to escape any characters, to avoid misinterpretation of special characters
while IFS=$'\t' read -r -a fields
do
    # Join relevant fields to create a grouping key
    # Fields 1(SampleType), 3(Time), and 4(Treatment) are collected in grouping_key
    grouping_key="${fields[1]}_${fields[3]}_${fields[4]}"
    # When comparing strings in bash, use double brackets
    if [[ ${grouping_key} != ${current_group} ]]
    then
        # Start a new group
        current_group=${grouping_key}
        echo "Group: ${current_group}" >> ${group_file}
    fi

    # Add the sample name to the current group, the sample name is in the 5th file. %% gets rid of the suffix.
    echo "SampleName: ${fields[5]%%_1.fq.gz}" >> ${group_file}
done < "${sorted_file}"

# Run the SAME while loop for a second time so that the issue with skipping "Group" rows in Step 2 can be overcome.
while IFS=$'\t' read -r -a fields
do
    # Join relevant fields to create a grouping key
    # Fields 1(SampleType), 3(Time), and 4(Treatment) are collected in grouping_key
    grouping_key="${fields[1]}_${fields[3]}_${fields[4]}"
    # When comparing strings in bash, use double brackets
    if [[ ${grouping_key} != ${current_group} ]]
    then
        # Start a new group
        current_group=${grouping_key}
        echo "Group: ${current_group}" >> ${group_file}
	echo "Group: ${current_group}" >> ${final_group_file}
    fi

    # Add the sample name to the current group
    echo "SampleName: ${fields[5]%%_1.fq.gz}" >> ${group_file}
    echo "SampleName: ${fields[5]%%_1.fq.gz}" >> ${final_group_file}
done < "${sorted_file}"

echo -e "The final grouping result is shown in ${final_group_file}"
echo -e "Note: to resolve an issue with skipping rows in Step 2, the collection of groups are delibrately ran TWICE. This is saved in ${group_file}, for data processing purpose only."

##### PROCESS Step 2 #####
# The output file of bedtool multicov will show the gene count in the final few columns, so we will just need to add those fields up to get average in STEP 5 (more on this later!)
bedfile="TriTrypDB-46_TcongolenseIL3000_2019.bed"

# Input file containing group information
group_file="groups_Tco2.txt"

# Loop through the groups in the input file
# IFS is set to an empty string, which means that no field splitting occurs and the entire line is treated as a single value
while IFS= read -r line
do
	# Check if the line starts with "Group:"
	if [[ $line == "Group:"* ]]
	then
		# Extract the group name by removing the prefix 'Group: '
		GroupName="${line#Group: }"
		echo "Group: ${GroupName}" # DEBUG LINE
		# Create an array to store sample names in the group
		sample_names=()
        # Read the sample names for this group, only executing the while loop if both conditions are true:
	# Condition 1: the line can be read
	# Condition 2: the line is not 'Group:'. This makes sure that it is the SampleName that it is reading 
        while IFS= read -r line && [[ $line != "Group:"* ]] 
	do
		if [[ $line == "SampleName:"* ]] 
		then
			# Extract the sample name and add it to the array, getting rid of 'SampleName: '
			sample_name="${line#SampleName: }"
			# Add them to the array, note the bracket used
			sample_names+=("$sample_name")
		fi
        done # End of sample_name collection

	# DEBUG LINE
	echo "Sample names: ${sample_names[@]}"

        # Run bedtools multicov for the samples in this group
	# Initialise an empty array 'bams' to store BAM file paths
	bams=()
	# '@' ensures that the entire array is used for the for loop (so we can put all of them as input for -bams 
        for sample_name in "${sample_names[@]}"
	do
		# Note the use of brackets for array
		bams+=("sorted_output_${sample_name}.bam")
		echo "End of BAM file paths collection, they are: ${bams}"
        done # End of BAM file paths collection

	# DEBUG LINE
	echo "BAM file paths: ${bams[@]}"

	# Remove the extra space at the end of this string variable, one % used to remove the shortest matching suffix (i.e.the last space)
        bams=${bams%" "}

        # Run bedtools multicov
	echo "Running bedtools multicov for group ${GroupName}" # DEBUG LINE
        bedtools multicov -bams "${bams[@]}" -bed "${bedfile}" > "${GroupName}_genecounts.cov"
        echo "Multicov completed for group ${GroupName}"
	fi
done < "${group_file}"

echo "STEP 4 completed, please proceed to STEP 5."
