#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
check='on'

# help message
help_message="
Wrapper to annotated peak BED files with Homer

usage:
    bash $(basename $0) [-options] --peaks <BED> --genome <STRING>
required arguments:
    -p|--peaks : BED file of peaks 
    -g|--genome : genome version [hg19,hg38,mm9,mm10]
optional arguments:
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -l|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on,off] (default = on)
    --depend : list of dependencies (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -p|--peaks)
            peaks=$2
            shift
            ;;
        -g|--genome)
            genome=$2
            shift
            ;;
        -o|--outdir)
            outdir="./$2"
            shift
            ;;
        -l|--logdir)
            logdir="./$2"
            shift
            ;;
        -n|--name)
            name=$2
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
            printf "ERROR: Unrecognised argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z ${peaks:-} ]]; then
    printf "\nERROR: --peaks argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${genome:-} ]]; then
    printf "\nERROR: --genome argument required\n"
    echo "$help_message"; exit 1
fi

# check files 
if [[ "${check:-}" = on ]]; then
    if [[ ! -r ${peaks} ]]; then
        printf "\nERROR: peaks file cannot be read: %s\n" $peaks
        echo "$help_message"; exit 1
    fi
fi

# set/check output directories
if [[ -z "${outdir:-}" ]]; then
    outdir="."
fi
mkdir -p $outdir
if [[ -z "${logdir:-}" ]]; then
    logdir="${outdir}"
fi
mkdir -p $logdir

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    echo "WARNING: name not set, extracting from: $peaks"
    name=$(basename $peaks)
    name=${name%%.*}
fi

# set log file names
logfile=${logdir}/${name}.$(basename "$0" .sh).log
scrfile=${logdir}/${name}.$(basename "$0" .sh).scr

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 20
		#SBATCH --mem=32gb
		#SBATCH --job-name=homer.$name
		#SBATCH -o ${logfile}
		${depend:-}

		source ~/miniconda3/etc/profile.d/conda.sh
		conda activate homer
		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		annotatePeaks.pl $peaks $genome -annStats ${outdir}/${name}_annStats.txt > \
			${outdir}/${name}_anno.txt
 
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
        
	EOS
) 
echo "$script" > $scrfile

# submit job
jobid=$(sbatch "$scrfile")

# echo job id and exit
echo "JOBID: $jobid"
exit 0 
