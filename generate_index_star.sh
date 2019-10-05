#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'

# help message
help_message="

purpose:
    Wrapper to run PBS script to generate STAR index
usage:
    bash $(basename "$0") [-options] -G <GTF>
required arguments:
    -G|--gtf : GTF file
    -F|--fasta : genome FASTA file
optional arguments:
    -n|--name : prefix for output files [string] (default = extracted from FASTA)
    -o|--outdir : output directory [string] (default = '.')
    -l|--logdir : output dir for log files [string] (default = --outdir)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list of PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # example depenency list: 'afterok:123456,afterok:123457'
output:
    # outputs STAR index files to --outdir & log file to --logdir

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -G|--gtf)
            gtf=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -l|--logdir)
            logdir=$2
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
		--depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arg
if [[ -z ${gtf:-} ]]; then
	printf "\nERROR: no --gtf argument provided\n" 
	echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
	printf "\nERROR: no --fasta argument provided\n"
	echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$gtf ]]; then
        printf "\nERROR: GTF file is not readable: %s/%s\n" $workdir $gtf
        echo "$help_message"; exit 1
    elif [[ ! -r $workdir/$fasta ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    fi
fi

# set name if not provided
if [[ -z ${name:-} ]]; then
    name=$(basename ${fasta})
    name=${name%%.*}
fi

# get basenames/prefixes
fasta_prefix=${fasta%.*}
fasta_base=$(basename $fasta)
gtf_base=$(basename $gtf)

# check if gzipped
if [[ ${fasta} == *".gz" ]]; then
    gzcmd="gzip -d ${fasta_base}"
    fasta_base=${fasta_base%.gz*}
fi

# create required dirs
mkdir -p ${workdir}/${logdir}
mkdir -p ${workdir}/${outdir}

# set log file names
scr_name=$(basename "$0" .sh)
scr_log=$workdir/$logdir/$name.$scr_name.scr.log
std_log=$workdir/$logdir/$name.$scr_name.std.log

# gen job script
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 20
		#SBATCH --mem=20G
		#SBATCH --job-name=index
		#SBATCH --output=$std_log

		# load modules
		module load STAR
	
		# copy files to scratch
		cp -L $workdir/$fasta_prefix* .
		cp -L $workdir/$gtf .

		# make outdir
		mkdir -p $name

		# run STAR
		STAR --runMode genomeGenerate \
			--runThreadN 20 \
			--genomeDir $name \
			--genomeFastaFiles $fasta_base \
			--sjdbGTFfile $gtf_base

		# copy output to outdir
		cp -r $name $workdir/$outdir
		
		ls -lhAR
	EOS
)
echo "$script" > $scr_log

# submit job
jobid=$(sbatch "$scr_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
