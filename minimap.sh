#!/bin/bash

### slurm script for minimap2

#SBATCH --job-name=mgCST_classifier                       # Job name
#SBATCH --output=                                         # Standard output log
#SBATCH --error=                                          # Standard error log
#SBATCH --ntasks=1                                        # Number of tasks
#SBATCH --mem=25G                                         # Memory per node
#SBATCH --partition=all                                   # Partition name

fastq_file=/path/concatenated_fastq_file.fastq
sam_file=experiment_name.sam
ref_gen=/path/ref_genome.fasta

minimap2 -ax map-ont -t 2 “$ref_gen” “$fastq_file” > “$sam_file”
