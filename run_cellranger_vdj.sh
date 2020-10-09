#!/user/bin/env bash
set -o pipefail
set -o nounset

# help message
help_message="
Wrapper to run deeptools bamCompare to generate tracks of difference in signal between two bam.

usage:
    bash $(basename "$0") [-options] -b <BAM>
required arguments:
    -f|--fastqs : Directory containing fastq from cellranger mkfastq
    -s|--sample : Sample ID to analyse
    -i|--id : Output sample ID 
    -r|--ref : Cellranger reference seq dir
optional arguments:
    -o|--outdir : output directory (default = pwd)
    -l|--logdir : output directory for log files (default = --outdir)
additional info:
    # all paths should be relative to the current working directory

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -f|--fastqs)
            fastqdir=$2
            shift
            ;;
        -s|--sample)
            sm_name=$2
            shift
            ;;
        -i|--id)
            name=$2
            shift
            ;;
        -r|--ref)
            ref=$2
            shift
            ;;
        -o|--outdir)
            outdir="./"
            shift
            ;;
        -l|--logdir)
            logdir="./$2"
            shift
            ;;
        *)
            echo "\nError: Unrecognised argument: %s %s" $1 $2
            echo "$help_message"; exit 1
        ;;
    esac
    shift
done

# set/check output directories
if [[ -z ${outdir:-} ]]; then
    outdir="."
fi
mkdir -p $outdir
if [[ -z ${logdir:-} ]]; then
    logdir=${outdir}
fi

# check required arguments
if [[ -z ${ref:-} ]]; then
    printf "\nERROR: no reference dir provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${sm_name:-} ]]; then
    printf "\nERROR: no sample name given\n"
    echo "$help_message"; exit 1
elif [[ -z ${name:-} ]]; then
    echo "\nERROR: no output ID given\n"
    echo "$help_message"; exit 1
elif [[ -z ${fastqdir} ]]; then
    echo "\nERROR: no fastq path given\n"
    echo "${help_message}"
fi

# set log names
scr_name=$(basename $0 .sh)
scr_log=${logdir}/${name}.${scr_name}.scr
std_log=${logdir}/${name}.${scr_name}.log

# run job
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 40
		#SBATCH --mem=60G
		#SBATCH --job-name=vdj_${name}
		#SBATCH --output=$std_log
	
		# load modules
		module load cellranger

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# generate coverage track
		cellranger vdj \
			--fastqs ${fastqdir} \
			--sample ${sm_name} \
			--id ${name} \
			--reference ${ref} \
			--localmem 60 \
			--localcores 40

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
	EOS
)
echo "$script" > $scr_log

# submit job
jobid=$(qsub "$scr_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
