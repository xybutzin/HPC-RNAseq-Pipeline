#!/bin/bash
# This script runs featureCounts on all sorted BAM files to qunatify reads for different annotation types (CDS, rRNA, gene, pseudogene).
# USAGE: sbatch fCount_on_allfiles.sh (if submitting as slurm job)
# USAGE: sh fCount_on_allfiles.sh (if run as a bash script)


#SBATCH --job-name=fCounts
#SBATCH --output=fCounts_%j.out # Standard output 
#SBATCH --error=fCounts_%j.err # Standard error log
#SBATCH --ntasks=1  # Run on a single CPU
#SBATCH --cpus-per-task=4 # Number of CPU cores per task (6)
#SBATCH --mem=32GB  # Job memory request (64)
#SBATCH --partition=compute
#SBATCH --time=2:00:00  # 24:00:00

# Paths to input files and directories
gff=/mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/genome/reference_genome/ecoli-CFP.gff3
bam_files=/mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/results/bowtie2/*.sortedByCoord.bam
fcounts_dir=/mmfs1/scratch/jacks.local/nicholas.butzin/RNAseq-Pipeline/results/featureCounts

# Get the current date (YYYYMMDD) for timestamping output files
current_date=$(date +%Y%m%d)

# Define output file paths for featureCounts results based on annotation types
featureCounts_cds_out=${fcounts_dir}/fcounts_CDS_${current_date}.txt
featureCounts_gene_out=${fcounts_dir}/fcounts_gene_${current_date}.txt
featureCounts_rRNA_out=${fcounts_dir}/fcounts_rRNA_${current_date}.txt
featureCounts_pseudogene_out=${fcounts_dir}/fcounts_pseudogene_${current_date}.txt

# Run featureCounts for different annotation types
# `-p`: pair-end reads
# `--countReadPairs`: Count read pairs instead of individual reads
# `-T 4`: Use 4 threads
# `-t <feature>`: Type of feature to count (e.g., CDS, rRNA)
# `-g gene`: Group counts by the gene attribute in the GFF file
# `-O`: Allow reads to be assigned to overlapping features
# `--fraction`: Assign a fractional count to each feature if a read overlaps multiple features
# `-a ${gff}`: Input GFF annotation file
# `-o <output_file>`: Path to the output file
# `<input_files>`: List of BAM files to process
featureCounts -p --countReadPairs -T 4 -t CDS -g gene -O --fraction -a ${gff} -o ${featureCounts_cds_out} ${bam_files}
featureCounts -p --countReadPairs -T 4 -t rRNA -g gene -O --fraction -a ${gff} -o ${featureCounts_rRNA_out} ${bam_files}
featureCounts -p --countReadPairs -T 4 -t gene -g gene -O --fraction -a ${gff} -o ${featureCounts_gene_out} ${bam_files}
featureCounts -p --countReadPairs -T 4 -t pseudogene -g gene -O --fraction -a ${gff} -o ${featureCounts_pseudogene_out} ${bam_files}

