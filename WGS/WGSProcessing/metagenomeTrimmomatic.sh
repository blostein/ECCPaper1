#! /bin/bash
# metagenomeTrimmomatic.sh
# William L. Close
# Pat Schloss Lab
# University of Michigan
#revised by F Blostein on May 13th 2020 to fit the meta_CAVITIES project parameters (namely, no head crop)
#revised so that zcat takes -f flag (allows uncompressed files) 

# Purpose: Quality trim and filter reads using Trimmomatic. Also removes 15 bases from start of both reads to account for adaptase (ideally would only remove from R2
# but Trimmomatic doesn't allow asymmetric trimming of one read).
# Usage: bash metagenomeTrimmomatic.sh FASTQ1 FASTQ2 OUTPUTDIR

##################
# Set Script Env #
##################

# Variables defined by user
FASTQ1=${1:?ERROR: Need to define FASTQ1.} # Forward read fastq
FASTQ2=${2:?ERROR: Need to define FASTQ2.} # Reverse read fastq

# Other variables
HEADCROP=${3:-15} # Use 15 as default value based on recommendations for Swift Accel NGS 1S Plus kit (trims extra bases from adaptase step)
ILLUMINACLIPADAPTERS=${4:-$(realpath $(which trimmomatic) | sed 's:trimmomatic$::')/adapters/TruSeq3-PE.fa} # Parses conda env PATH to use TruSeq3-PE adapters unless changed by the user
#OUTDIR=/scratch/bfoxman_root/bfoxman/blostein/meta_CAVITIES/Files/CleanPipe/PostTrim
OUTDIR=${5:?ERROR: Please provide output directory.} #output directory



#########################################################
# Removing Adapters and Trimming Reads with Trimmomatic #
#########################################################

echo PROGRESS: Trimming raw sequence reads.

# Create output folder if it doesn't already exist
echo ${OUTDIR}
#mkdir -p "${OUTDIR}"/

# Pulling FASTQ[12] file names minus the extensions for labelling Trimmomatic output
BASENAME=$(echo "${FASTQ1}" | sed 's/.*\/\(.*\)_R._001\.fastq\.gz/\1/')
#BASENAME=$(echo "${FASTQ1}" | sed 's/\(.*\)-.*/\1/')
FILENAME1=$(echo "${FASTQ1}" | sed 's/.*\/\(.*\)_001\.fastq\.gz/\1/')
#FILENAME1=$(echo "${FASTQ1}" | sed 's/\(.*\)\..*/\1/')
FILENAME2=$(echo "${FASTQ2}" | sed 's/.*\/\(.*\)_001\.fastq\.gz/\1/')
#FILENAME2=$(echo "${FASTQ2}" | sed 's/\(.*\)\..*/\1/')

# Running Trimmomatic
TrimmomaticPE \
	"${FASTQ1}" \
	"${FASTQ2}" \
    "${OUTDIR}"/"${FILENAME1}"_paired.fastq.gz \
    "${OUTDIR}"/"${FILENAME1}"_unpaired.fastq.gz \
    "${OUTDIR}"/"${FILENAME2}"_paired.fastq.gz \
    "${OUTDIR}"/"${FILENAME2}"_unpaired.fastq.gz \
    ILLUMINACLIP:"${ILLUMINACLIPADAPTERS}":2:30:10:4:true \
    SLIDINGWINDOW:4:20 MINLEN:75 HEADCROP:"${HEADCROP}"

# Combining unpaired reads into single output file
mv "${OUTDIR}"/"${FILENAME1}"_unpaired.fastq.gz "${OUTDIR}"/"${BASENAME}"_unpaired.fastq.gz && gunzip "${OUTDIR}"/"${BASENAME}"_unpaired.fastq.gz
zcat -f "${OUTDIR}"/"${FILENAME2}"_unpaired.fastq.gz >> "${OUTDIR}"/"${BASENAME}"_unpaired.fastq && gzip "${OUTDIR}"/"${BASENAME}"_unpaired.fastq && rm "${OUTDIR}"/"${FILENAME2}"_unpaired.fastq.gz

#fastqc
mkdir -p "${OUTDIR}"/FASTQC

fastqc "${OUTDIR}"/"${FILENAME1}"_paired.fastq --outdir="${OUTDIR}"/FASTQC
fastqc "${OUTDIR}"/"${FILENAME2}"_paired.fastq --outdir="${OUTDIR}"/FASTQC
fastqc "${OUTDIR}"/"${BASENAME}"_unpaired.fastq.gz --outdir="${OUTDIR}"/FASTQC

