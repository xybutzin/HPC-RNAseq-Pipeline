#!/bin/bash
# This script processes paired-end FASTQ files (R1, R2) through a RNA-seq pipeline.
# It performs quality control (FastQC), alignment (bowtie2), and feature counting (samtools, and featureCounts).
# Usage: sh rnaseq_analysis_on_input_file.sh <name of R1 fastq file> <name of R2 fastq file>

# Change directories to /mmfs1/scratch/jacks.local/xuan.butzin where all analysis results will be stored
cd /mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline
echo "Current working directory: $(pwd)"

# Initialize variables to store the names of the input fastq file
fq1=$1  # Input FASTQ file for R1 (forward reads)
fq2=$2  # Input FASTQ file for R2 (reverse reads)

# Extract the base name from filename for naming outputs
samplename=`basename ${fq1} _R1.fastq.gz`
echo "Sample name is ${samplename}"

# Specify the number of CPU cores to use
cores=6

# Define the paths to the genome index files and gene annotation file
genome=/mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/genome/index_genome/listed_genome/ecoli_CFP
gff=/mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/genome/reference_genome/ecoli-CFP.gff3

# Set up output filenames and locations
fastqc_out=results/fastqc/
align_out=results/bowtie2/
align_out_sam=$align_out${samplename}_aligned.sam
align_out_bam=$align_out${samplename}_aligned.bam
align_out_sorted_bam=$align_out${samplename}_aligned.sortedByCoord.bam
fcounts_out=results/featureCounts

# Create all the output directories if they don't exist (-p)
mkdir -p $fastqc_out $align_out $fcounts_out

# Load requred software modules (use version numbers)
module load fastqc/0.12.1 
#module load gcc/11.2.0
module load bowtie2/2.5.2
module load samtools/1.19

# Run FastQC for quality control on the input FASTQ files and save results to the designated directory
echo "Starting QC for ${samplename}"
fastqc ${fq1} ${fq2} -o ${fastqc_out}

# Perform alignment with Bowtie2, specifying the genome index, input FASTQ files, and output SAM file
echo "Starting bowtie2 for ${samplename}"
bowtie2 -p 4 --local --time -x ${genome} -1 ${fq1} -2 ${fq2} -S ${align_out_sam}

# Convert the SAM file to a BAM file using Samtools
echo "Starting samtools to convert sam to bam for ${samplename}"
samtools view -hbS ${align_out_sam} > ${align_out_bam}

# Sort the BAM file by genomic coordinates using Samtools
echo "Sorting bam file for ${samplename}"
samtools sort ${align_out_bam} -o ${align_out_sorted_bam}

# Index the sorted BAM file using Samtools
echo "Indexing for ${samplename}"
samtools index ${align_out_sorted_bam}

echo "Done"

