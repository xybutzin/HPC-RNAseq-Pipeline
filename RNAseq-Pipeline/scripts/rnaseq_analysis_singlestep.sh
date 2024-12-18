#!/bin/bash

# This script automates RNA-seq analysis by iterating over paired-end FASTQ files(R1 and R2)
# in a specified directory, performing quality control, alignment, sorting and feature count all in one step

# TODO: right now the session needs to stay open in order for the featureCount check loop to continuously check. 
# The script will stop if the session is closed
 
# Usage: sh rnaseq_analysis_singlestep.sh
# Ensure the directory path and file patterns match your data structure


# Set up output filenames and locations
output_dir="/mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/results"
fastqc_out=$output_dir/fastqc/
align_out=$output_dir/bowtie2/
fcounts_out=$output_dir/featureCounts


# Define the paths to the genome index files and gene annotation file
genome=/mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/genome/index_genome/listed_genome/ecoli_CFP
gff=/mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/genome/reference_genome/ecoli-CFP.gff3

# Get the list of R1 FASTQ files and count the total number of samples
fq1_files=(/mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/data/*_R1.fastq.gz)
total_samples=${#fq1_files[@]}
echo "Found $total_samples total samples"


# Loop through all R1 FASTQ files in the specified directory
for fq1 in "${fq1_files[@]}"; do
    # Derive the corresponding fq2 file name by replacing "_R1.fastq.gz" with "_R2.fastq.gz"
    fq2="${fq1/_R1.fastq.gz/_R2.fastq.gz}"
    
    echo "Submitting job for: $fq1 and $fq2"
    # Submit a job to SLURM for RNA-seq analysis using sbatch
    sbatch -p compute -t 12:00:00 -c 4 --job-name rnaseq-workflow --mem 64G \
      --wrap="sh /mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/scripts/rnaseq_analysis_on_input_file.sh ${fq1} ${fq2}"
    
    # Add a brief delay to avoid overloading the job scheduler
    sleep 2
done


# Continuous check to ensure all the BAM files are ready before starting featureCounts
echo "Checking if all sorted BAM files are ready for featureCounts..."

# Loop to check every 15 minutes (900 seconds) until all BAM files are ready
while true; do
    # Count the number of sorted BAM files present
    current_bam_count=$(find "$align_out" -type f -name "*_aligned.sortedByCoord.bam" | wc -l)

    if [ "$current_bam_count" -eq "$total_samples" ]; then
        echo "All $current_bam_count BAM files are ready. Proceeding to featureCounts."
        break
    else
        echo "Only $current_bam_count out of $total_samples BAM files are ready. Rechecking in 20 seconds..."
        sleep 20
    fi
done


sbatch fCount_on_allfiles.sh

