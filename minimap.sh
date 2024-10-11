#!/bin/bash

### slurm script for minimap2

#SBATCH --job-name=mgCST_classifier                       # Job name
#SBATCH --output=/local/scratch/amaros/logs/slurm/R.out   # Standard output log
#SBATCH --error=/local/scratch/amaros/logs/slurm/R.err    # Standard error log
#SBATCH --ntasks=1                                        # Number of tasks
#SBATCH --mem=25G                                         # Memory per node
#SBATCH --partition=all                                   # Partition name

fastq_file=path_to_fastq_directory/concatenated_fastq_file.fastq
sam_file=experiment_name.sam
gen_ref=path_to_ref_genome.fasta

minimap2 -ax map-ont -t 2 “$ref_file” “$fastq_file” > “$sam_file”
