#!/bin/bash

#set input directory
#set input and output directory
in=/scratch/bfoxman_root/bfoxman/blostein/meta_CAVITIES/Files/CleanPipe/PostTrim/

# get count of files in this directory
NUMFILES=$(ls -1 ${in}*R1*.fastq.gz | wc -l)

# subtract 1 as we have to use zero-based indexing (first element is 0)
ZBNUMFILES=$(($NUMFILES - 1))

# submit array of jobs to SLURM
if [ $ZBNUMFILES -ge 0 ]; then
  sbatch --array=1-58%5 array_decontam.pbs
else
  echo "No jobs to submit, since no input files in this directory."
fi
