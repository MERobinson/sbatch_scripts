#!/usr/bin/env bash
set -o pipefail
set -o nounset 

# default arg
check='on'

# help message
help_message="
Wrapper to call differential footprints from DNase data with Wellington

usage:
    bash $(basename "$0") [-options] -r <BED> -b1 <BAM> -b2 <BAM>
required arguments:
    -r|--regions : accessible regions to profile [BED]
    -t|--treatment : aligned reads for treatment/test [BAM]
    -c|--control : aligned reads for control [BAM] 
optional arguments:
    -o|--outdir : output directory name (default = PWD)
    -n|--name : name of comparison & outdir subfolder (default = <tname>_vs_<cname>)
    -tn|--tname : treatment name (default = extracted from treatment)
    -cn|--cname : control name (defailt = extracted from control)
    --fdr : FDR cutoff (default = 0.01)
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
		-t|--treatment)
		    treatment=$2
		    shift
		    ;;
        -c|--control)
            control=$2
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
		-tn|--tname)
		    tname=$2
		    shift
		    ;;
        -cn|--cname)
            cname=$2
            shift
            ;;
        --fdr)
            fdr="-fdr $2"
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
		    printf "ERROR: Unrecognised argument: %s %s\n" $1 $2
		    echo "$help_message"; exit 1
		    ;;
	esac
	shift
done

# check required args
if [[ -z ${regions:-} ]]; then
    printf "\nERROR: --region argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${treatment:-} ]]; then
    printf "\nERROR: --treatment argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${control:-} ]]; then
    printf "\nERROR: --control argument required\n"
    echo "$help_message"; exit 1
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

# set optional args
if [[ -z ${tname:-} ]]; then
    tname=$(basename "${treatment}")
    tname=${tname%%.*}
fi
if [[ -z ${cname:-} ]]; then
    cname=$(basename ${control})
    cname=${cname%%.*}
fi
if [[ -z ${name:-} ]]; then
    name="${tname}_vs_${cname}"
fi

# check files
if [[ ${check} = 'on' ]]; then
    if [[ ! -r ${regions} ]]; then
        printf "\nERROR: input file cannot be read: %s\n" $regions
        echo "$help_message"; exit 1
    elif [[ ! -r ${treatment} ]]; then
        printf "\nERROR: bam file cannot be read: %s\n" $treatment
        echo "$help_message"; exit 1
    elif [[ ! -r ${control} ]]; then
        printf "\nERROR: bam file cannot be read: %s\n" $control
        echo "$help_message"; exit 1
    fi
fi

# set commands
bed_cmd=("zcat -f ${regions} | "
         "awk 'BEGIN {OFS = \"\t\"}"
         "{if (\$3 - \$2 <  100) \$2 = \$2 - 50; \$3 = \$3 + 50; print \$1, \$2, \$3}'"
         "> ${name}.preprocessed.bed")
pydnase_cmd=("wellington_bootstrap.py -A ${fdr_limit:-} -p 12" 
             "${teatment} ${control} preprocessed.bed"
             "${name}/${tname}_specific_fp.txt"
             "${name}/${cname}_specific_fp.txt")

# set log file names
logfile=${logdir}/${name}.$(basename "$0" .sh).log
scrfile=${logdir}/${name}.$(basename "$0" .sh).scr

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 12
		#SBATCH -mem=24gb
		#SBATCH --job-name=diffFP.${name}
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
echo "$script" > $scrfile
sbatch "$scrfile"
exit 0
