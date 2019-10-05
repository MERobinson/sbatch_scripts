#!/user/bin/env bash
set -o pipefail
set -o nounset

# default arg
format='bedgraph'
outdir=''
extend=100
quality=10
dedup=yes
binsize=100
libtype='SE'
check='yes'
scr_name=$(basename "$0")

# help message
help_message="
usage:
    bash $scr_name [-options] -b <BAM>
purpose:
    # Wrapper to run deeptools bamCoverage to generate genome tracks.
required arguments:
    -b|--bam : Input bam file
    -F|--fasta : Reference FASTA file
optional arguments:
    -o|--outdir : output directory for track files (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output file (default = BAM filename)
    -f|--format : output format [bedgraph,bigwig] (default = bedgraph)
    -t|--libtype : library type - paired or single [PE|SE] (default = SE)
    -c|--center : whether to center reads [INT] (default = NULL)
    -g|--genome : genome version - used in trackline if provided (default = NULL)
    -e|--extend : bp extension for generating coverage track [INT] (default = 100)
    -q|--quality : min alignment quality to filter reads (default = 10)
                   (to disable simply set to 0)
    -d|--dedup : whether to filter duplicates [yes,no] (default = yes)
    -s|--binsize : binning for track output (default = 100)
    --smooth : window size for smoothing (default = NULL)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL) 
additional info:
    # all paths should be relative to the current working directory
    # FASTA file is used to obtain a list of canonical chromosomes for filtering
      the BAM prior to generating tracks (regex for canonical = chr[0-9MXY]+)
    # Note if --libtype PE is set, extension and center sizes is ignored and
      actual isize from paired mappings used 

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam=$2
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
        -t|--libtype)
            libtype=$2
            shift
            ;;
        -c|--center)
            center=$2
            shift
            ;;
        --smooth)
            smooth=$2
            shift
            ;;
        -f|--format)
            format=$2
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
        --check)
            check=$2
            shift
            ;;
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        *)
            printf "\nError: Illegal argument: %s %s" $1 $2
            echo "$help_message"; exit 1
        ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arguments
if [[ -z ${bam:-} ]]; then
    printf "\nERROR: no BAM file provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: no FASTA file provided\n"
    echo "$help_message"; exit 1
fi

# set base
bam_base=$(basename "$bam")
bam_pref=${bam_base%%.*}

# if no name provided extract from BAM
if [[ -z "${name:-}" ]]; then
    name=${bam_pref}
fi

# set tmpdir
tmpdir="tmp_${name}_trackgen"
if [[ -e $tmpdir ]]; then
    printf "\nERROR: tmpdir already exists: %s\n" $tmpdir
    exit 1
fi
mkdir $tmpdir

# set optional argument
if [[ ${dedup} = yes ]]; then
    dedup_arg="-F 1024"
fi
if [[ ${libtype} = 'PE' ]]; then
    ext_arg="--extendReads"
else
    ext_arg="--extendReads ${extend}"
fi
if [[ -n ${center:-} ]]; then
    center_arg="--centerReads"
fi
if [[ -n ${smooth:-} ]]; then
    smooth_arg="--smoothLength $smooth"
fi

# set format dependent arg
if [[ $format = 'bedgraph' ]]; then
	extension="bedGraph"
    trackline="track type=$extension name=$name visibility=2 windowingFunction=mean autoScale=off"
    trackline="${trackline} smoothingWindow=10 viewLimits=0:75 maxHeightPixels=0:75:150"
    if [[ -n ${db:-} ]]; then trackline="${trackline} db=$db"; fi
    format_arg="cat $name.$extension | grep -E \"^chr[0-9XYM]+\b\" > tmp.txt"
    format_arg="${format_arg:-}; echo $trackline | cat - tmp.txt > $name.$extension"
elif [[ $format = 'bigwig' ]]; then
    extension='bw'
else
    printf "\nERROR: format not recognised, must be either bedgraph or bigwig\n"
    echo "$help_message"; exit 1
fi

# get list of canonical chromosomes
chr_list=$(cat $fasta | grep -Eo "^>chr[0-9MXY]+\b" | \
           grep -Eo "chr[0-9XY]+"| tr "\n" " ") 

# create required dirs
mkdir -p $outdir
mkdir -p $logdir

# set log file names
scr_name=${scr_name%.*}
std_log=$logdir/$name.$scr_name.std.log
scr_log=$logdir/$name.$scr_name.scr.log

# compile commands
bamcov_cmd=("bamCoverage ${ext_arg:-} ${center_arg:-} ${smooth_arg:-}"
                "--normalizeUsing RPKM --ignoreForNormalization chrX chrM"
                "--numberOfProcessors 12 --outFileFormat ${format}"
                "--bam ${tmpdir}/input.filt.bam --binSize ${binsize}"
                "--outFileName ${tmpdir}/${name}.${extension}")

# run job
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 12
		#SBATCH --mem=18G
		#SBATCH --job-name=trackgen
		#SBATCH --output=$std_log

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# load modules
		module load samtools
		source activate deeptools

		# copy resource files to scratch
		mkdir -p $tmpdir
		cp ${bam}* $tmpdir/ 

		# filter bam prior to making track
		samtools view -b ${dedup_arg:-} -q $quality \
			$tmpdir/$bam_base $chr_list > $tmpdir/input.filt.bam
		samtools index $tmpdir/input.filt.bam

		# generate coverage track
		${bamcov_cmd[@]} 

		${format_arg:-}	
		mv $tmpdir/$name.$extension $outdir/
		ls -lhAR $tmpdir
		rm -r $tmpdir
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
	EOS
)
echo "$script" > $scr_log

# submit job
jobid=$(qsub "$scr_log")

# echo jobid and exit
echo "JOBID: $jobid"
exit 0
