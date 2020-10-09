#!/usr/bin/env bash
set -o pipefail
set -o nounset 

# default arg
outdir=''
check='on'
libtype='A'
mode='quasi'

# help message
help_message="
Wrapper to quantify transcripts with Salmon v0.12.0

usage:
    bash $(basename "$0") [-options] -r1 <FASTQ|BAM> -i <INDEX|FASTA> 
required arguments:
    -r1|--reads1 : input file(s) - either FASTQ or BAM dependening on mode
    -i|--index : reference - either salmon index or FASTA depending on mode
optional arguments:
    -r2|--reads2 : 2nd FASTQ file (if PE & quasi mode) 
    -m|--mode : quasi or alignment based [quasi|align] (default = quasi)
    -l|--libtype : library types [string - see salmon man for options] (default = A)
    -n|--name : output filename prefix (default = extracted from input file)
    -o|--outdir : output directory name (default = PWD)
    --logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on,off] (default = on)
additional info:
    # if quasi mode --reads should be FASTQ files and --index salmon index
    # if align mode --reads should be aligned BAM and --index transcripts FASTA
    # if multiple runs/reps - read files should be given as comma sep list

"

# parse command line arguments
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -r1|--reads1)
            r1=$2
            shift
            ;;
        -i|--index)
            index=$2
            shift
            ;;
        -r2|--reads2)
            r2=$2
            shift
            ;;
        -m|--mode)
            mode=$2
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
        -l|--libtype)
            libtype=$2
            shift
            ;;
        --logdir)
            logdir="./$2"
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        *)
            printf "\nERROR: Unrecognised argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required args
if [[ -z ${r1:-} ]]; then
    printf "\nERROR: --reads1 argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${index:-} ]]; then
    printf "\nERROR: --index argument required\n"
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

# split potential multi runs
IFS=',' read -r -a r1_array <<< ${r1}
if [[ -n ${r2:-} ]]; then
    IFS=',' read -r -a r2_array <<< ${r2:-}
fi

# set optional args
if [[ -z ${name:-} ]]; then
    name=$(basename "${r1}")
    name=${name%%.*}
fi
if [[ $mode = 'quasi' ]]; then
    r1_arg="-r ${r1_array[@]}"
	idx_arg="-i ${index}"
elif [[ $mode = 'align' ]]; then
    r1_arg="-a ${r1_array[@]}"
	idx_arg="-t ${index}"
else
    printf "\nERROR: --mode argument not recognised: %s\n" $mode 
    echo "$help_message"; exit 1
fi
if [[ -n ${r2_array:-} ]]; then
    r1_arg="-1 ${r1_array[@]}"
    r2_arg="-2 ${r2_array[@]}"
fi

# check files
if [[ ${check} = 'on' ]]; then
    for r in ${r1_array[@]}; do  
        if [[ ! -r ${r} ]]; then 
            printf "\nERROR: input file cannot be read: %s/%s\n" $r
            echo "$help_message"; exit 1
        fi
    done
    if [[ ! -r ${index} ]]; then
        printf "\nERROR: index cannot be read: %s/%s\n" $index
        echo "$help_message"; exit 1
    fi
    if [[ -n ${r2:-} ]]; then
        for r in ${r2_array[@]}; do
            if [[ ! -r ${r} ]]; then 
                printf "\nERROR: input file cannot be read: %s/%s\n" $r
                echo "$help_message"; exit 1
            fi
        done
    fi
fi

# set commands
salmon_cmd=("salmon quant ${idx_arg} -l ${libtype}" 
            "${r1_arg} ${r2_arg:-} -o ${outdir}/${name}"
            "--validateMappings -p 12")

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$logdir/$name.$scr_name.log
scr_log=$logdir/$name.$scr_name.scr

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#SBATCH --time=12:00:00
		#SBATCH -N 1
		#SBATCH -n 12
		#SBATCH --mem=24gb
		#SBATCH --job-name=$name.salmon
		#SBATCH --output=$std_log

		# load modules
		source ~/miniconda3/etc/profile.d/conda.sh
		conda activate salmon

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# run Salmon
		${salmon_cmd[@]}
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
	EOS
)
echo "$script" > $scr_log

# submit job, echo id and exit
jobid=$(sbatch "$scr_log")
echo "JOBID: $jobid"
exit 0
