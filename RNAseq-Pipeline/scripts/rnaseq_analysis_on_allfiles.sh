#! /bin/bash

# This script automates RNA-seq analysis by iterating over paired-end FASTQ files(R1 and R2)
# in a specified directory and submitting jobs to a SLURM cluster
 
# Usage: sh rnaseq_analysis_on_allfiles.sh
# Ensure the directory path and file patterns match your data structure

# Loop through all R1 FASTQ files in the specified directory
for fq1 in /mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/data/C*_R1.fastq.gz
do
    # Derive the corresponding fq2 file name by replacing "_R1_001.fastq.gz" with "_R2_001.fastq.gz"
    fq2="${fq1/_R1.fastq.gz/_R2.fastq.gz}"
    
    echo "Submitting job for: $(basename "$fq1") and $(basename "$fq2")"
    # Submit a job to SLURM for RNA-seq analysis using sbatch
    sbatch -p compute -t 12:00:00 -c 4 --job-name rnaseq-workflow --mem 64G --wrap="sh /mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/scripts/rnaseq_analysis_on_input_file.sh ${fq1} ${fq2}"
    
    # Add a brief delay to avoid overloading the job scheduler
    sleep 2
done