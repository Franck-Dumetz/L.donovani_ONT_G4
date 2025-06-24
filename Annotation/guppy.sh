#!/bin/bash
#$ -cwd
#$ -P dserre-lab
#$ -o qsubOutFinal.txt
#$ -e qsubErrFinal.txt
#$ -l mem_free=5G
#$ -sync y
#$ -t 1-337
#$ -tc 10
#$ -V

### Basecalling script for ONT direct RNA sequencing using Guppy

export LD_LIBRARY_PATH=/usr/local/packages/python-3.8/lib:$LD_LIBRARY_PATH

fast5_dir=path_to_fast5_files
output_dir=path_to_output_directory

/usr/local/packages/guppy-6.4.2_gpu/bin/guppy_basecaller --input_path "$fast5_dir" --save_path "$output_dir" --config rna_r9.4.1_70bps_hac.cfg --min_qscore 7 --records_per_fastq 10000000 --gpu_runners_per_device 8 --recursive -x 'cuda:all' --num_callers 8
