#!/usr/bin/env bsh
set -o errexit
set -o pipefail
set -o nounset

# default arg
check='on'
chimeric='off'
twopass='off'

# help message
help_message="
Wrapper to align reads with STAR and produce BAM & counts.

usage:
    bash $(basename "$0") [-options] -f1 <FASTQ1> -I <INDEX>
required arguments:
    -fq1|--fastq1 : input FASTQ file
    -I|--index : STAR index directory
optional arguments:
    -fq2|--fastq2 : mate FASTQ file if PE [FASTQ] (default = NULL)
    -n|--name : prefix for output files [string] (default = extracted from fastq)
    -o|--outdir : output directory [string] (default = '.')
    -l|--logdir : output dir for log files [string] (default = --outdir)
    --chimeric : whether to detect chimeric reads [on|off] (default = 'off')
    --twopass : whether to run in 2pass mode [on|off] (default = 'off')
    --sjdb : comma sep list of SJ.out.tab files to pass to STARs
             --sjdbFileChrStartEnd argument (default = NULL)
    --clip3 : n of bases to clip from 3' end (default = NULL)
    --check : whether to check input files [on|off] (default = 'on')
    --depend : list of PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # example depenency list: 'afterok:123456,afterok:123467'

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -fq1|--fastq1)
            fq1=$2
            shift
            ;;
		-fq2|--fastq2)
            fq2=$2
            shift
            ;;
        -I|--index)
            index=$2
            shift
            ;;
        -n|--name)
            name=$2
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
        --chimeric)
            chimeric=$2
            shift
            ;;
        --twopass)
            twopass=$2
            shift
            ;;
        --sjdb)
            sjdb=$2
            shift
            ;;
        --clip3)
            clip3=$2
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

# set/check output directories
if [[ -z ${outdir:-} ]]; then
    outdir="."
fi
mkdir -p $outdir
if [[ -z ${logdir:-} ]]; then
    logdir=${outdir}
fi
mkdir -p $logdir
tmpdir=$(mktemp -d)

# check required arg
if [[ -z ${fq1:-} ]]; then
    printf "\nERROR: no --fastq1 argument provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${index:-} ]]; then
    printf "\nERROR: no --index argument provided\n"
    echo "$help_message"; exit 1
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
    name=$(basename "$fq1")
    name=${name%%.*}
fi

# set compression arg
fq_ext=${fq1##*.}
if [[ $fq_ext = "gz" ]]; then
    compress_arg="--readFilesCommand zcat"
    fq_ext="fastq.gz"
else
    fq_ext="fastq"
fi

# check if multiple fastq per read
IFS=',' read -r -a fq1_array <<< $fq1
if [[ ${#fq1_array[@]} -gt 1 ]]; then
    echo "Merging input R1 fastq"
    check='off'
    fq1_merge="cat ${fq1_array[@]##*/} > ${tmpdir}/${name}.R1_merge.${fq_ext}"
    fq1="${tmpdir}/${name}.R1_merge.${fq_ext}"
fi
if [[ -n ${fq2:-} ]]; then
    echo "Mate input - running in PE mode"
    IFS=',' read -r -a fq2_array <<< $fq2
    if [[ ${#fq2_array[@]} -gt 1 ]]; then
        echo "Merging input R2 fastq"
        fq2_merge="cat ${fq2_array[@]##*/} > ${tmpdir}/${name}.R2_merge.${fq_ext}"
        fq2="${tmpdir}/${name}.R2_merge.${fq_ext}"
    fi
else
    echo "No mate input - running in SE mode"
fi

# check files unless flagged
if [[ $check = 'on' ]]; then
    if [[ ! -r $fq1 ]]; then
        printf "\nERROR: FASTQ file is not readable: %s\n" $fq1
        echo "$help_message"; exit 1
    elif [[ -n ${fq2:-} ]] && [[ ! -r ${fq2:-} ]]; then
        printf "\nERROR: FASTQ file is not readable: %s\n" ${fq2:-}
        echo "$help_message"; exit 1
    elif [[ ! -r $index ]]; then
        printf "\nERROR: Index folder is not readable: %s\n" $index
        echo "$help_message"; exit 1
    fi
fi

# check optional flags
if [[ $chimeric = 'on' ]]; then
    chimeric_arg="--chimSegmentMin 15"
fi
if [[ $twopass = 'on' ]]; then
    twopass_arg="--twopassMode Basic"
fi
if [[ -n ${sjdb:-} ]]; then
    IFS=',' read -r -a sjdb_array <<< "$sjdb"
    sjdb_arg="--sjdbFileChrStartEnd ${sjdb_array[@]##*/}"
fi
if [[ -n ${clip3:-} ]]; then
    clip3_arg="--clip3pNbases ${clip3}"
fi

# set commands
star_command=("STAR --runMode alignReads"
                "--runThreadN 20"
                "--genomeDir ${index}"
                "--readFilesIn ${fq1} ${fq2:-}"
                "--outFilterType BySJout"
                "--outFilterMultimapNmax 20"
                "--alignSJoverhangMin 8"
                "--alignSJDBoverhangMin 1"
                "--outFilterMismatchNoverLmax 0.04"
                "--alignIntronMin 20"
                "--alignIntronMax 1000000"
                "--alignMatesGapMax 1000000"
                "--outSAMstrandField intronMotif"
                "--quantMode TranscriptomeSAM GeneCounts"
                "--readNameSeparator _"
                "--outFileNamePrefix ${outdir}/${name}.star."
                "--outBAMsortingThreadN 20"
                "--outSAMtype BAM SortedByCoordinate"
                "${compress_arg:-}"
                "${chimeric_arg:-}"
                "${twopass_arg:-}"
                "${sjdb_arg:-}"
                "${clip3_arg:-}")

# set log file names
scr_name=$(basename "$0" .sh)
std_log=${logdir}/${name}.${scr_name}.log
scr_log=${logdir}/${name}.${scr_name}.scr

# run PBS script
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -n 20
		#SBATCH -N 1
		#SBATCH --mem=60G
		#SBATCH --job-name=$name.star
		#SBATCH --output=$std_log

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# load modules
		source ~/miniconda3/etc/profile.d/conda.sh
        conda activate star
		
		# merge fastq if multiple
		${fq1_merge:-} 
		${fq2_merge:-}

		# run STAR
		${star_command[@]} 
		
		# index bam
		mv ${outdir}/${name}.star.Aligned.sortedByCoord.out.bam \
			${outdir}/${name}.star.bam 
		samtools index ${outdir}/${name}.star.bam
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` 
		EOS
) 
echo "$script" > $scr_log

# submit job
jobid=$(sbatch "$scr_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
