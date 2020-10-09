#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
check='on'
chrlist='hs'

# help message
help_message="

Takes an aligned BAM then filters & shifts tags.

usage:
    bash $(basename $0 .sh) [-options] -i <BAM>
required arguments:
    -i|--input : input BAM file [BAM]
    -bl|--blacklist : input blacklist regions [BED]
optional arguments:
    -chr|--chrlist : chr to include - either species or comma sep list of chr
                     [mm|hs|list] (default = 'hs')
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on|off] (default = on)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified
    # if --chrlist = mm/hs, chr included = chr[0-9XY]+

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--input)
            input=$2
            shift
            ;;
        -bl|--blacklist)
            blacklist=$2
            shift
            ;;
        -chr|--chrlist)
            chrlist=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -q|--qcdir)
            qcdir=$2
            shift
            ;;
        -l|--logdir)
            logdir=$2
            shift
            ;;
        -n|--name)
            name=$2
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
            printf "ERROR: Unrecognised argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z ${input:-} ]]; then
    printf "\nERROR: --input argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${blacklist:-} ]]; then
    printf "\nERROR: --blacklist argument required\n"
    echo "$help_message"; exit 1
fi

# check files 
if [[ "${check:-}" = on ]]; then
    if [[ ! -r ${input} ]]; then
        printf "\nERROR: Input BAM cannot be read: %s/%s\n" $input
        echo "$help_message"; exit 1
    elif [[ ! -r ${blacklist} ]]; then
        printf "\nERROR: Input blacklist BED file cannot be read\n" $blacklist
        echo "$help_message"; exit 1
    fi
fi

# set output directories
if [[ -z ${outdir:-} ]]; then
    outdir='.'
fi
mkdir -p $outdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
mkdir -p $logdir
if [[ -z "${qcdir:-}" ]]; then
    qcdir=$outdir
fi
mkdir -p $qcdir
tmpdir=$(mktemp -d)

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    name=$(basename ${input})
    name=${name%%.*}
fi

# set chr list
if [[ ${chrlist} = 'hs' ]]; then
    chr_array=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" 
               "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" 
               "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
elif [[ ${chrlist} = 'mm' ]]; then
    chr_array=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9"
               "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19"
               "chrX" "chrY")
else
    IFS=',' read -r -a chr_array <<< ${chrlist}
fi

# set commands
chrfilt_cmd=("samtools view -q 10 -o ${tmpdir}/${name}.chrfilt.bam"
             "${input} ${chr_array[@]}")
blfilt_cmd=("bedtools intersect -v -a ${tmpdir}/${name}.chrfilt.bam"
            "-b ${blacklist} > ${tmpdir}/${name}.blfilt.bam")
filt_cmd=("samtools view -F 1548 -u ${tmpdir}/${name}.blfilt.bam |"
           "samtools sort -@ 8 -n -o ${tmpdir}/${name}.dedup.bam -")
fixmate_cmd=("samtools fixmate -r ${tmpdir}/${name}.dedup.bam ${tmpdir}/${name}.fix.bam")
bedpe_cmd=("bedtools bamtobed -bedpe -mate1 -i ${tmpdir}/${name}.fix.bam |"
           "gzip -nc > ${tmpdir}/${name}.bedpe.gz")
shift_cmd=("zcat -f ${tmpdir}/${name}.bedpe.gz |"
           "awk 'BEGIN {OFS = \"\t\"}"
           "{if (\$9 == \"+\") {\$2 = \$2 + 4}"
           "else if (\$9 == \"-\") {\$3 = \$3 - 5}"
           "if (\$10 == \"+\") {\$5 = \$5 + 4}"
           "else if (\$10 == \"-\") {\$6 = \$6 - 5} print \$0}' |"
           "gzip -nc > ${outdir}/${name}.tn5.bedpe.gz")
bedtota_cmd=("zcat ${outdir}/${name}.tn5.bedpe.gz |"
             "awk 'BEGIN{OFS=\"\t\"}"
             "{printf \"%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n\","
             "\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' |"
             "gzip -nc > ${outdir}/${name}.tn5.tagAlign.gz")

# set log file names
logfile=$logdir/$name.$(basename ${0} .sh).log
scrfile=$logdir/$name.$(basename ${0} .sh).scr

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=48:00:00
		#SBATCH -N 1
		#SBATCH -n 1
		#SBATCH --mem=18G
		#SBATCH --job-name=${name}.atac
		#SBATCH --output=${logfile}

		# load modules
		source ~/miniconda3/etc/profile.d/conda.sh
		conda activate sambedtools
		mkdir -p ${tmpdir}
		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
		
		# filter
		echo "Filtering BAM"
		${chrfilt_cmd[@]}
		${blfilt_cmd[@]}
		${filt_cmd[@]}

		# fix mate
		echo "Fixing mates"
		${fixmate_cmd[@]}

		# bedpe conversion
		echo "Converting to BEDPE"
		${bedpe_cmd[@]} 

		# tag shifting
		echo "Tn5 shifting BEDPE reads"
		${shift_cmd[@]}

		# tag conversion
		echo "Converting filtered, shifted BEDPE to tagAlign"
		${bedtota_cmd[@]} 

		# sort and index bam
		samtools sort -o ${outdir}/${name}.filt.bam ${tmpdir}/${name}.fix.bam
		samtools index ${outdir}/${name}.filt.bam
		samtools flagstat ${outdir}/${name}.filt.bam > ${qcdir}/${name}.filt.flagstats.txt

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
	EOS
) 
echo "$script" > ${scrfile}

# submit job
jobid=$(sbatch "$scrfile")
echo "$jobid"
exit 0 
