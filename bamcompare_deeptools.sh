#!/user/bin/env bash
set -o pipefail
set -o nounset

# default arg
format='bigwig'
outdir=${PWD}
extend=200
quality=10
dedup='yes'
binsize=200
mode='log2'

# help message
help_message="
Wrapper to run deeptools bamCompare to generate tracks of difference in signal between two bam.

usage:
    bash $(basename "$0") [-options] -b <BAM>
required arguments:
    -t|--test : BAM file for test sample
    -c|--cntrl : BAM file for control sample
    -F|--fasta : Reference FASTA file
optional arguments:
    -o|--outdir : output directory for track files (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output file (default = BAM filename)
    -f|--format : output format [bedgraph,bigwig] (default = bedgraph)
    -m|--mode : type of comparison [log2|ratio|subtract|mean] (default = 'log2')
    -g|--genome : genome version - used in trackline if provided (default = NULL)
    -e|--extend : bp extension for generating coverage track (default = 200)
    -q|--quality : min alignment quality to filter reads (default = 10)
                   (to disable simply set to 0)
    -d|--dedup : whether to filter duplicates [yes,no] (default = yes)
    -s|--binsize : binning for track output (default = 200)
additional info:
    # all paths should be relative to the current working directory
    # FASTA file is used to obtain a list of canonical chromosomes for filtering
      the BAM prior to generating tracks (regex for canonical = chr[0-9MXY]+)

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -t|--test)
            test_bam=$2
            shift
            ;;
        -c|--cntrl)
            control_bam=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
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
        -f|--format)
            format=$2
            shift
            ;;
        -m|--mode)
            mode=$2
            shift
            ;;
        -g|--genome)
            db=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -e|--extend)
            extend=$2
            shift
            ;;
        -q|--quality)
            quality=$2
            shift
            ;;
        -d|--dedup)
            dedup=$2
            shift
            ;;
        -s|--binsize)
            binsize=$2
            shift
            ;;
        --smooth)
            smooth=$2
            shift
            ;;
        *)
            echo "\nError: Illegal argument: %s %s" $1 $2
            echo "$help_message"; exit 1
        ;;
    esac
    shift
done

# set dir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
mkdir -p ${outdir}
mkdir -p ${logdir}

# check required arguments
if [[ -z ${test_bam:-} ]]; then
    printf "\nERROR: no test BAM file provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: no FASTA file provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${control_bam:-} ]]; then
    printf "\nERROR: no control BAM file provided\n"
    echo "$help_message"; exit 1
fi

# set base
test_base=$(basename "$test_bam")
control_base=$(basename "$control_bam")

# if no name provided extract from BAM
if [[ -z ${name:-} ]]; then
    name=${test_base%%.*}
fi

# set tmpdir
tmpdir="tmp_${name}_bamcompare"
if [[ -e ${tmpdir} ]]; then
    printf "\nERROR: tmpdir already exists: %s\n" $tmpdir
    exit 1
fi
mkdir $tmpdir

# set dedup argument
if [[ $dedup = yes ]]; then
    dedup_arg="-F 1024"
fi
if [[ -n ${smooth:-} ]]; then
    smooth_arg="--smoothLength $smooth"
fi

# set format dependent arg
if [[ $format = 'bedgraph' ]]; then
	extension="bedGraph"
elif [[ $format = 'bigwig' ]]; then
    extension='bw'
else
    printf "\nERROR: format not recognised, must be either bedgraph or bigwig\n"
    echo "$help_message"; exit 1
fi
    
# get list of canonical chromosomes
chr_list=$(cat $fasta | grep -Eo "^>chr[0-9MXY]+\b" | \
           grep -Eo "chr[0-9XY]+"| tr "\n" " ") 

# set log names
scr_name=$(basename $0 .sh)
scr_log=$logdir/$name.$scr_name.scr.log
std_log=$logdir/$name.$scr_name.std.log

# run job
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 20
		#SBATCH --mem=20G
		#SBATCH --job-name=bamcomp
		#SBATCH --output=$std_log
	
		# load modules
		module load samtools
		source activate deeptools    

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# copy resource files to scratch
		cp ${test_bam}* ${tmpdir}/
		cp ${control_bam}* ${tmpdir}/

		# filter bam prior to making track
		samtools view -b ${dedup_arg:-} -q ${quality} \
			${tmpdir}/${test_base} ${chr_list} > ${tmpdir}/${name}.test.bam
		samtools view -b ${dedup_arg:-} -q ${quality} \
			${tmpdir}/${control_base} ${chr_list} > ${tmpdir}/${name}.control.bam
		samtools index ${tmpdir}/${name}.test.bam
		samtools index ${tmpdir}/${name}.control.bam

		# generate coverage track
		bamCompare \
			--bamfile1 ${tmpdir}/${name}.test.bam \
			--bamfile2 ${tmpdir}/${name}.control.bam \
			--operation ${mode} \
			--extendReads ${extend} \
			--binSize ${binsize} \
			--scaleFactorsMethod readCount \
			--ignoreForNormalization chrX \
			--numberOfProcessors 20 \
			--outFileFormat ${format} \
			--outFileName ${name}.${extension} ${smooth_arg}

		# move output to outdir
		mv ${name}.${extension} $outdir/

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
		ls -lhAR ${tmpdir}
		rm -r ${tmpdir}
	EOS
)
echo "$script" > $scr_log

# submit job
jobid=$(qsub "$scr_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
