#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
check='yes'
date=`date '+%Y-%m-%d %H:%M:%S'`

# help message
help_message="
Merge BAM files from separate runs and generate metrics for merged alignments.

usage:
    bash $(basename $0) [-options] -b <BAM1,BAM2> -f <FASTA>
required arguments:
    -b|--bam_list : comma separated list of input BAM files to merge
    -f|--fasta : whole genome FASTA
optional arguments:
    -n|--name : name prefix for output files (default = extracted from first bam)
    -o|--outdir : outut directory for merged BAM (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --depend : comma-sep list of job dependencies - format afterok:<jobid>,afterok:<jobid>
    --check : whether to check files prior to running job [yes,no] (default = yes)
additional_info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # bam/log/qc output directories inherit from --outdir unless specified

"

# parse command line arg
while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
        -b|--bam_list)
            bam_list=$2
            shift
            ;;
        -f|--fasta)
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
        -q|--qcdir)
            qcdir=$2
            shift
            ;;
        --depend)
            depend="#SBATCH --dependency=$2"
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required arg provided
if [[ -z "${bam_list:-}" ]]; then
    printf "\nERROR: --bam_list argument required\n"
    echo "$help_message"; exit 1
else
    IFS="," read -r -a bam_array <<< "$bam_list"
    for bam in ${bam_array[@]}; do
        if [[ $check = "yes" ]] && [[ ! -r $bam ]]; then
            printf "\nERROR: BAM file not readable: %s\n" $bam
            echo "$help_message"; exit 1
        else
            merge_input="${merge_input:-}I=$bam "
        fi
    done
fi
if [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
else
    if [[ $check = 'yes' ]] && [[ ! -r $fasta ]]; then 
        printf "\nERROR: BAM file not readable: %s\n" $bam
        echo "$help_message"; exit 1
    fi
fi

# if no name provided extract from first bam
if [[ -z ${name:-} ]]; then
    name="${bam_array[0]%%.*}.merged"
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
if [[ -z "${qcdir:-}" ]]; then
    qcdir="${outdir}"
fi
mkdir -p $qcdir
tmpdir=$(mktemp -d)

# set commands
merge_cmd=("picard -Xmx16g MergeSamFiles"
            "${merge_input}"
            "OUTPUT=${outdir}/${name}.bam"
            "TMP_DIR=${tmpdir}")
metrics_cmd=("picard -Xmx16g CollectAlignmentSummaryMetrics"
             "R=${fasta}"
             "I=${outdir}/${name}.bam"
             "O=${qcdir}/${name}.alignment_summary_metrics.txt"
             "TMP_DIR=${tmpdir}")

# set log file names
logfile=${logdir}/${name}.$(basename "$0" .sh).log
scrfile=${logdir}/${name}.$(basename "$0" .sh).scr

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 12
		#SBATCH --mem=18gb
		#SBATCH --job-name=${name}.merge
		#SBATCH -o ${logfile}
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# load modules
		source ~/miniconda3/etc/profile.d/conda.sh
		conda activate sambedpicard
		mkdir -p ${tmpdir}

		# merge 
		${merge_cmd[@]}

		# alignment metrics
		${metrics_cmd[@]}

		# copy final bam to outdir
		samtools index ${outdir}/${name}.bam

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
	EOS
)
echo "$script" > $scrfile
sbatch "$scrfile"
exit 0
