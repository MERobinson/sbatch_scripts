#!/usr/bin/env bash
set -o pipefail
set -o nounset 

# default arg
check='on'

# help message
help_message="
Wrapper to call footprints from DNase data with Wellington

usage:
    bash $(basename "$0") [-options] -r <BED> -b <BAM>
required arguments:
    -r|--regions : accessible regions to profile [BED]
    -b|--bam : BAM files to find footprints in [BAM]
optional arguments:
    -n|--name : output filename prefix (default = extracted from input file)
    -o|--outdir : output directory name (default = PWD)
    --fdrlimit : minimum pval to consider sig for FDR (default = -20)
    --logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on,off] (default = on)
    --depend : list of PBS dependencies (default = NULL)
additional info:

"

# parse command line arguments
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-r|--regions)
		    regions=$2
		    shift
		    ;;
		-b|--bam)
		    bam=$2
		    shift
		    ;;
		-o|--outdir)
		    outdir=$2
		    shift
		    ;;
		-n|--name)
		    name=$2
		    shift
		    ;;
        --fdrlimit)
            fdrlimit="-fdrlimit $2"
            shift
            ;;
        --logdir)
            logdir=$2
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
		    echo "ERROR: Unrecognised argument: %s %s" $1 $2
		    echo "$help_message"; exit 1
		    ;;
	esac
	shift
done

# check required args
if [[ -z ${regions:-} ]]; then
    printf "\nERROR: --region argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${bam:-} ]]; then
    printf "\nERROR: --bam argument required\n"
    echo "$help_message"; exit 1
fi

# set optional args
if [[ -z ${name:-} ]]; then
    name=$(basename "${bam}")
    name=${name%%.*}
fi
if [[ -z "${outdir:-}" ]]; then
    outdir="${name}"
fi
mkdir -p $outdir
if [[ -z "${logdir:-}" ]]; then
    logdir="${outdir}"
fi
mkdir -p $logdir

# check files
if [[ ${check} = 'on' ]]; then
    if [[ ! -r "${regions}" ]]; then
        printf "\nERROR: input file cannot be read: %s\n" $regions
        echo "$help_message"; exit 1
    elif [[ ! -r "${bam}" ]]; then
        printf "\nERROR: index cannot be read: %s\n" $bam
        echo "$help_message"; exit 1
    fi
fi

# set commands
bed_cmd=("zcat -f ${regions} | "
         "awk 'BEGIN {OFS = \"\t\"}"
         "{if (\$3 - \$2 <  100) \$2 = \$2 - 50; \$3 = \$3 + 50; print \$1, \$2, \$3}'"
         "> ${name}.preprocessed.bed")
pydnase_cmd=("wellington_footprints.py -A ${fdr_limit:-} -p 12" 
             "-o ${name} ${name}.preprocessed.bed ${bam} ${outdir}/${name}")

# set log file names
logfile=${logdir}/${name}.$(basename "$0" .sh).log
scrfile=${logdir}/${name}.$(basename "$0" .sh).scr

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#SBATCH --time=48:00:00
		#SBATCH -N 1
		#SBATCH -n 12
		#SBATCH --mem=32gb
		#SBATCH --job-name=$name.fp
		#SBATCH -o ${logfile}
		${depend:-}

		# load modules
		source ~/miniconda3/etc/profile.d/conda.sh
		conda activate pydnase
		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# pre process bed
		${bed_cmd[@]}

		# run wellington
		mkdir ${outdir}/${name}
		${pydnase_cmd[@]}

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
	EOS
)
echo "$script" > ${scrfile}

# submit job, echo id and exit
sbatch "$scrfile"
exit 0
