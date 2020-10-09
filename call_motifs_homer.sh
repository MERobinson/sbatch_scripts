#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default args
check='on'
size='-size 200'
denovo='off'
scan='off'

# help message
help_message="
Wrapper to run HOMERs findMotifsGenome.pl

usage:
    bash $(basename $0) [-options] -b <BED> -g <GENOME>
required args:
    -b|--bed : BED file of regions to analyse 
    -g|--genome : genome version [STRING]
optional args:
    --size : frag size [INT|'given'] (default = 200)
    -bg|--background : background regions [BED] (default = NULL)
    --mknown : known motif file to use in place of in-built DB (default = NULL)
    --denovo : whether to run denovo motif discovery [on|off] (default = off)
    --scan : only scan for known motifs NOT run enrichment [on|off] (default = off)
    -n|--name : output filename prefix (default = FASTA prefix)
    -o|--outdir : output directory for results (default = PWD)
    -l|--logdir : output directory for logs (default = --outdir)
    --check : whther to check input files [on|off] (default = on)
    --depend : list of dependencies (default = NULL)
additional info:
    # all paths should be relative to working dir
    # check and depend options used for job scheduling
    # log output dir inherits from --outdir unless specified
    # see HOMER docs for more info
    # scan mode requires known motif file to be preovided 

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bed)
            bed="$2"
            shift
            ;;
        -g|--genome)
            genome="$2"
            shift
            ;;
        -bg|--background)
            bg="$2"
            shift
            ;;
        --size)
            size="-size ${2}"
            shift
            ;;
        --mknown)
            mknown="$2"
            shift
            ;;
        --denovo)
            denovo="$2"
            shift
            ;;
        --scan)
            scan="$2"
            shift
            ;;
        -o|--outdir)
            outdir="$2"
            shift
            ;;
        -l|--logdir)
            logdir="$2"
            shift
            ;;
        -n|--name)
            name="$2"
            shift
            ;;
        --check)
            check="$2"
            shift
            ;;
        --depend)
            depend="#SBATCH --dependency=$2"
            shift
            ;;
    esac
    shift
done

# check required args
if [[ -z ${bed:-} ]]; then
    printf "\nERROR: --bed argument is required\n"
    echo "$help_message"; exit 1
elif [[ -z ${genome:-} ]]; then
    printf "\nERROR: --genome argument if required\n"
    echo "$help_message"; exit 1
fi

# check files
if [[ $check = 'on' ]]; then
    if [[ ! -r $bed ]]; then
        printf "\nERROR: input BED file cannot be read: %s\n" $bed
        echo "$help_message"; exit 1
    elif [[ -n ${bg:-} ]] && [[ ! -r $bg ]]; then
        printf "\nERROR: input motif file cannot be read: %s\n" $background
        echo "$help_message"; exit 1
    fi
fi

# set/check output directories
if [[ -z "${outdir:-}" ]]; then
    outdir="."
fi
mkdir -p "${outdir}"
if [[ -z "${logdir:-}" ]]; then
    logdir="${outdir}"
fi
mkdir -p "${logdir}"

# set name
if [[ -z ${name:-} ]]; then
    name=$(basename "$bed")
    name=${name%%.*}
fi

# handle optional args
if [[ -n ${bg:-} ]]; then
    bg_arg="-bg ${bg}"
fi
if [[ -n ${mknown:-} ]]; then
    mknown_arg="-mknown ${mknown}"
fi
if [[ ${denovo} = 'on' ]]; then
    denovo_arg="-nomotif"
fi
if [[ ${scan} = 'on' ]]; then
    if [[ -z ${mknown:-} ]]; then
        printf "\n\nERROR: scanning mode (-find) requires known motif to be provided\n"
        exit 1
    fi
    scan_arg="-find ${mknown} > ${outdir}/${name}.motif_occurence.txt"
fi

# set log file names
logfile=${logdir}/${name}.$(basename "$0" .sh).log
scrfile=${logdir}/${name}.$(basename "$0" .sh).scr

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 12
		#SBATCH --mem=24gb
		#SBATCH --job-name=homer.$name
		#SBATCH -o ${logfile}
		${depend:-}

		source ~/miniconda3/etc/profile.d/conda.sh
		conda activate homer
		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		findMotifsGenome.pl ${bed} ${genome} ${outdir} ${size}\
			${bg_arg:-} ${mknown_arg:-} ${denovo_arg:-} ${scan_arg:-}
 
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`    
	EOS
)
echo "$script" > ${scrfile}

# submit job
jobid=$(sbatch "$scrfile")
echo "${jobid}"
exit 0
