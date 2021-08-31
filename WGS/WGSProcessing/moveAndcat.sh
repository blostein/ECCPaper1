#! /bin/bash
#moveAndcat.sh
#F Blostein
#University of Michigan

#Purpose: Move trimmed, cleaned and decontaminated files from individual sample folders to another folder and concatonate them
#for use with humann2
#Usage: bash moveAndcat.sh

##########################
#Set Script Env#
##########################

#Variables 
my_dir=/scratch/bfoxman_root/bfoxman/blostein/meta_CAVITIES/Files/CleanPipe/PostHR/
new_dir=/scratch/bfoxman_root/bfoxman/blostein/meta_CAVITIES/Humann2/test_files

cd ${my_dir}
for basename in $(ls *.fq.gz | cut -f1 -d_ | sort -u)
do
    cat "$basename"*R1*.fq.gz "$basename"*R2*.fq.gz > "${basename}_cat.fq.gz"
    mv "${basename}_cat.fq.gz" ${new_dir}
done
