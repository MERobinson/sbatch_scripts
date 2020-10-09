#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
sradir="$HOME/sradata/sra"

# help message
help_message="
Wrapper to pre-fetch and extract FASTQ data using SRA-tools.

usage:
    bash $(basename $0) [-options] -s <SRA_accession>
required arguments:
    -s|--sra : SRA accession number
optional arguments:
    -n|--name : name prefix for output files (default = SRA accession)
    -o|--outdir : output directory (default = PWD)
    -ld|--logdir : output directory for log files (default = --outdir)
    -sd|--sradir : directory setup for sra toolkit (default = ~/sradata/sra)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    # all paths (except sradir) should be given relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified
    # make sure sradir points to the directory setup for SRA toolkit

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -s|--sra)
            sra=$2
            shift
            ;;
        -o|--outdir)
            outdir="./$2"
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -ld|--logdir)
            logdir="./$2"
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        --depend)
            depend="#SBATCH --dependency=$2"
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check/parse required arguments
if [[ -z ${sra:-} ]]; then
    printf "\nERROR: --sra argument required.\n"
    echo "$help_message"; exit 1
fi

# set/check output directories
if [[ -z ${outdir:-} ]]; then
    outdir="."
fi
mkdir -p $outdir
if [[ -z ${logdir:-} ]]; then
    logdir=${outdir}
fi
mkdir -p $logdir

# set name if not provided
if [[ -z "${name:-}" ]]; then
    name=${sra}
fi

# set log file names
scrfile=${logdir}/${name}.$(basename ${0} .sh).scr
logfile=${logdir}/${name}.$(basename ${0} .sh).log

# job
script=$(cat <<- EOS
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 1
		#SBATCH --mem=6gb
		#SBATCH --job-name=${name}.sra
		#SBATCH -o ${logfile}
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
	
		# load modules
		module load SRA-Toolkit

		# prefetch
		if [[ ! -e ${sradir}/${sra}.sra ]]; then
			echo "Pre-fetching SRA:" 
			prefetch --max-size 500G ${sra}
		else
			echo "SRA already exists"
		fi

		# dump
		echo "Extracting FASTQ:"
		fastq-dump -I -W -B --skip-technical --outdir ${outdir} \
			--origfmt --gzip --split-files ${sradir}/${sra}.sra
		rename ${sra} ${name} ${outdir}/${sra}*.fastq.gz

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
	EOS
)
echo "${script}" > ${scrfile}

# submit job
jobid=$(sbatch "${scrfile}")
echo "${jobid}"
exit 0
