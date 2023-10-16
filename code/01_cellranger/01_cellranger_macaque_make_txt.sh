#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=2G
#$ -N BLA_mktxt_macaque
#$ -o logs/Macaque/mk_text.txt
#$ -e logs/Macaque/mk_text.txt
#$ -m e


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"



# Specify the directory where the sample folders are located
TARGET_DIR="/dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/raw-data/FASTQ/Macaque" 
output_file="01_cellranger_macaque.txt"

# Clearing the output file if it exists
echo -n > "$output_file"


# Rename folders based on the prefix
for dir in "$TARGET_DIR"/sample*; do
    # Check if it's a directory
    if [ -d "$dir" ]; then
        # Extract the prefix name from one of the files in the directory
        for file in "$dir"/*fastq.gz; do
            # Using 'awk' to get the desired prefix
            prefix=$(echo "$file" | awk -F'/' '{print $NF}' | awk -F'_S[0-9]*_L[0-9]*_[IR][0-9]*_001.fastq.gz' '{print $1}')
            
            # Rename the directory to the prefix
            mv "$dir" "$TARGET_DIR/$prefix"
            
            # Break after the first file as we need only one prefix per directory
            break
        done
    fi
done

# Now, collect the names of the renamed folders and write them to the output file
for new_dir in "$TARGET_DIR"/*; do
    if [ -d "$new_dir" ]; then
        echo $(basename "$new_dir") >> "$output_file"
    fi
done

echo "**** Job ends ****"
date