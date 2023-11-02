##### STEP 6 #####
echo -e "Finally, we are in STEP 6! The final step of this ICA."
echo -e "Our goal here is to calculate the fold change for group-wise comparisons, and we will do it for all groups"
echo -e "i.e. group 1 vs group 2, group 1 vs group 3, group 2 vs group 3 etc etc"

##### PROCESS Step 1 #####
# Creating an array of file names corresponding to the average files generated from STEP5
average_file=(*"_average.txt")

# Our strategy here is to use an awk script within a for loop.
# Use a for loop, call each item in the average_file array each time, as "file1", and record the iteration key as i. i increases by 1 for each iteration.
for ((i=0; i < ${#average_file[@]}; i++))
do
        file1="${average_file[i]}"
        # Set the group name
        group_name1="${file1%%_average.txt}"
        # Setting up a stacked for loop so that group_name1 has something to compared against, in this case, we are comparing group_name1 VS group_name2
        # j is always 1 bigger than i, and this will come in handy later to avoid comparing group2 VS group1, which are duplicates to values in group1 VS group2
        for ((j=i+1; j < ${#average_file[@]}; j++))
        do
                file2="${average_file[j]}"
                group_name2="${file2%%_average.txt}"
                # We only want to compare group1 vs group 2, but not against itself, nor group2 vs group1
                if [[ "${i}" < "${j}" ]]
                then
                        echo -e "Calculating fold change for ${group_name1} VS ${group_name2}"
                        # Use awk to get us the average value for both groups
                        avg1=$(awk '{print $3}' "${file1}")
                        avg2=$(awk '{print $3}' "${file2}")

                        #Don't forget about the gene names and descriptions, only one gene_info variable is needed because everything is already indexed within ${average_file}
                        gene_info=$(awk '{print $1, $2}' "${file1}")

                        # We want to avoid dividing by 0
                        if [[ ${avg2} != 0 ]]
                        then
                                # scale=5 will return 5 decimal places, bc allows calculations, or else the result is just a fraction
                                fold_change=$(echo "scale=5; ${avg1} / ${avg2}" | bc)
                        else
                                # Handle division by zero if avg2 is zero
                                fold_change=0
                        fi
                        # Using paste function combined with echo -e. Paste command can be used to concatenate files line by line. This is to make sure that the columns are added as new columns, not simply appended to the first column as new rows.
                        # The output of this will still contain fold_change added both as a new column and a new row. I spent way too much time debugging this, but realised the files generated could still be used for the next step
                        paste <(echo -e "${gene_info}") <(echo -e "${fold_change}") -d $'\t' > "${group_name1}_VS_${group_name2}_fold.txt"
                fi
        done
done

# Filtering out those unwanted rows
fold_file=(*"_fold.txt")
for fold_file in "${fold_file[@]}"
do
        # This part filters the unwanted rows
        awk 'NF == 3' ${fold_file} > "${fold_file}_filtered"
        # Renaming them to overwrite, by force, or else the user will have to type y for many times!
        mv -f "${fold_file}_filtered" "${fold_file}"
        echo -e "Filtering complete for ${fold_file}"
        # This part sorts by the 3rd column in descending orders
        # -V so that it sorts by 'natural sort' but not by the first character
        sort -k3,3r -V ${fold_file} > "${fold_file}_sorted.txt"
        echo -e "Sorting complete for ${fold_file}"
done

echo -e "The final results are stored in <group_name1>_VS_<group_name2>_fold.txt_sorted.txt"
echo -e "As a reminder, the grouping results which contains the sequence details are stored in final_groups_Tco2.txt, NOT groups_Tco2.txt"
echo -e "To view the first few lines of your result, type:"
echo -e "head <group_name1>_VS_<group_name2>_fold.txt_sorted.txt"
echo -e "For example: Sample Type:Clone1; Time: 0 Hours; Treatment: Uninduced VSType:Clone1; Time: 24 Hours; Treatment: Uninduced "
echo -e "head Clone1_0_Uninduced_VS_Clone1_24_Uninduced_fold.txt_sorted.txt
 "
echo -e "------THIS IS THE END OF STEP 6 AND THE END OF THE ICA------"
