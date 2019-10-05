#!/usr/bin/env bash
set -o pipefail
set -o nounset 

# default arg
workdir=$PWD
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
    --depend : list of PBS dependencies (default = NULL)
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
            outdir=$2
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

# split potential multi runs
IFS=',' read -r -a r1_array <<< ${r1}
if [[ ${#r1_array[@]} -gt 1 ]]; then
    for r in ${r1_array[@]}; do
        r1_bn="${r1_bn:-} $(basename ${r})"
        r1_cp="${r1_cp:-} cp ${workdir}/${r}* .;"
    done
else
    r1_bn=$(basename ${r1})
    r1_cp="cp ${workdir}/${r1}* ."
fi
if [[ -n ${r2:-} ]]; then
    IFS=',' read -r -a r2_array <<< ${r2:-}
    for r in ${r2_array[@]}; do
        r2_bn="${r2_bn:-} $(basename ${r})"
        r2_cp="${r2_cp:-} cp ${workdir}/${r}* .;"
    done
else
    r2_bn=$(basename ${r2})
    r2_cp="cp ${workdir}/${r2}* ."
fi

# set optional args
if [[ -z ${name:-} ]]; then
    name=$(basename "${r1}")
    name=${name%%.*}
fi
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
if [[ $mode = 'quasi' ]]; then
    r1_arg="-r ${r1_bn}"
	idx_arg="-i $(basename ${index})"
elif [[ $mode = 'align' ]]; then
    r1_arg="-a ${r1_bn}"
	idx_arg="-t $(basename ${index})"
else
    printf "\nERROR: --mode argument not recognised: %s\n" $mode 
    echo "$help_message"; exit 1
fi
if [[ -n ${r2:-} ]]; then
    r1_arg="-1 ${r1_bn}"
    r2_arg="-2 ${r2_bn}"
fi
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# check files
if [[ ${check} = 'on' ]]; then
    for r in ${r1_array[@]}; do  
        if [[ ! -r ${workdir}/${r} ]]; then 
            printf "\nERROR: input file cannot be read: %s/%s\n" $workdir $r
            echo "$help_message"; exit 1
        fi
    done
    if [[ ! -r "${workdir}/${index}" ]]; then
        printf "\nERROR: index cannot be read: %s/%s\n" $workdir $index
        echo "$help_message"; exit 1
    fi
    if [[ -n ${r2:-} ]]; then
        for r in ${r2_array[@]}; do
            if [[ ! -r ${workdir}/${r} ]]; then 
                printf "\nERROR: input file cannot be read: %s/%s\n" $workdir $r
                echo "$help_message"; exit 1
            fi
        done
    fi
fi

# set commands
salmon_cmd=("salmon quant ${idx_arg} -l ${libtype}" 
            "${r1_arg} ${r2_arg:-} -o salmon_out/${name}"
            "--validateMappings -p 12")

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$workdir/$logdir/$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=10:00:00
		#PBS -l select=1:ncpus=12:mem=10gb
		#PBS -j oe
		#PBS -N salmon.$name
		#PBS -q med-bio
		#PBS -o ${std_log}
		${depend:-}

		# load modules
		module load anaconda3/personal
		source activate salmon-0.12.0

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy input & index to scratch
		${r1_cp}
		${r2_cp:-}
		cp -rL ${workdir}/${index}* .

		# run Salmon
		${salmon_cmd[@]} &>> $out_log

		# copy results to output directory
		cp -r salmon_out/* $workdir/$outdir/
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR &>> $out_log
		ls -lhAR 
		cp -r $out_log $workdir/$logdir/
	EOS
)
echo "$script" > $pbs_log

# submit job, echo id and exit
jobid=$(qsub "$pbs_log")
echo "JOBID: $jobid"
exit 0
