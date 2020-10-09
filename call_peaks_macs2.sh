#!/bin/bash
set -o pipefail
set -o nounset

# default arg
genome='hs'
check='on'
force='off'
bdg='off'
broad='off'
nomodel='off'
summits='off'
dedup='on'
chrfilt='on'
qual='10'
format='BAM'

# help message
help_message="
Wrapper to call peaks with MACS2 

Usage:
    bash $(basename "$0") [-options] -t <BAM>
required arguments:
    -t|--test : test sample(s) [BAM]
    -F|--fasta : genome seq [FASTA]
optional arguments:
    -c|--control : control sample(s) [BAM]
    -o|--outdir : output dir for macs files (default = '.')
    -l|--logdir : output dir for log files (default = --outdir)
    -q|--qcdir : output dir for qc files (default = ---outdir)
    -g|--genome : genome size [hs,mm,<numeric>] (default = hs)
    -n|--name : prefix for output files (default = test sample name)
    -e|--extsize : run macs2 with provided extension (default = unset)
    -f|--format : format of test file [BAM|BED|...] (default = 'BAM')
    --shift : shift size [INT] (default = NULL)
    --broad : whether to run in broad peak mode [on|off] (default = off)
    --bdg : whether to generate bedgraph tracks [on|off] (default = off)
    --nomodel : whether to run in nomodel mode [on|off] (default = off)
    --pval : pvalue cutoff [INT] (default NULL)
    --summits : whether to call summits [on|off] (default = off)
    --dedup : whether to filter dups [on|off] (default = on)
    --chrfilt : whether to filter non-standard chr [on|off] (default = on)
    --qual : mapping quality filter score [numeric] (default = 10)
    --blfilt : blacklist regions to filter out [BED] (default = none)
    --force : whether to re-run if output already exist [on|off] (default = off)
    --check : whether to check files prior to running [on|off] (default = on)
    --depend : list of dependencies for PBS script (e.g. afterok:012345,afterok:012346)
additional info:
    # all paths should be given relative to working directory
    # --outdir will be created, if not pre-existing
    # if --control is provided, will be used as background for peak calling,
      otherwise no control sample is used in macs
    # if not provided, name is extracted from test filename up to first period
    # test/control can be given as single inputs or matched comma separated lists

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -t|--test)
            test=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -c|--control)
            control=$2
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
        -q|--qcdir)
            qcdir="./$2"
            shift
            ;;
        -g|--genome)
            genome=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -e|--extsize)
            extsize=$2
		    shift
		    ;;
        -f|--format)
            format="$2"
            shift
            ;;
        --nomodel)
            nomodel=$2
            shift
            ;;
        --shift)
            shift_arg="--shift $2"
            shift
            ;;
        --tsize)
            tagsize="--tsize $2"
            shift
            ;;
        --pval)
            pval_arg="--pvalue $2"
            shift
            ;;
        --summits)
            summits=$2
            shift
            ;;
        --broad)
            broad=$2
            shift
            ;;
        --bdg)
            bdg=$2
            shift
            ;;
        --dedup)
            dedup=$2
            shift
            ;;
        --chrfilt)
            chrfilt=$2
            shift
            ;;
        --qual)
            qual=$2
            shift
            ;;
        --blfilt)
            blfilt=$2
            shift
            ;;
        --force)
            force=$2
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
            printf "\nERROR: Unrecognised argument: %s %s\n" $1 $2
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
if [[ -z "${qcdir:-}" ]]; then
    qcdir=${outdir}
fi
mkdir -p $qcdir
tmpdir=$(mktemp -d)

# check required arguments
if [[ -z "${test:-}" ]]; then
    printf "\nERROR: No test BAM provided\n"
    echo "$help_message"; exit 1
elif [[ -z "${fasta:-}" ]]; then
    printf "\nERROR: No FASTA file provided\n"
    echo "$help_message"; exit 1
elif [[ ${check} = 'on' ]] && [[ ! -r $fasta ]]; then
    printf "\nERROR: FASTA cannot be read: %s\n" $fasta
    echo "$help_message"; exit 1
fi

# check optional arg
if [[ ${bdg} = 'on' ]]; then
    bdg_arg="--bdg --SPMR"
fi
if [[ ${summits} = 'on' ]]; then
    summits_arg="--call-summits"
fi
if [[ ${broad} = 'on' ]]; then
    broad_arg="--broad"
fi
if [[ ${nomodel} = 'on' ]]; then
    model_arg="--nomodel"
fi
if [[ -n ${extsize:-} ]]; then
    ext_arg="--extsize ${extsize}"
    model_arg="--nomodel"
fi

# split input lists to arrays
IFS="," read -r -a test_array <<< $test
if [[ -n ${control:-} ]]; then
    IFS="," read -r -a cntrl_array <<< $control
else
    echo "WARNING: no control samples input"
fi

# set name if not provided
if [[ -z ${name:-} ]]; then
    name_array=(${test_array[@]/%.*/})
    name_array=(${name_array[@]/#*\//})
else
	IFS="," read -r -a name_array <<< $name
fi

# get canonical chr from fasta
chr_list=$(cat $fasta | grep -Eo "^>(chr)?[0-9XY]+\b" | \
           grep -Eo "(chr)?[0-9XY]+"| tr "\n" " ")
if [[ $chr_list =~ "21" ]] && [[ $genome = "mm" ]]; then
    echo "ERROR: genome set as mouse but fasta contains chr21"
    exit 1
elif [[ ! $chr_list =~ "21" ]] && [[ $genome = "hs" ]]; then
    echo "ERROR: genome set as human but fasta doesn't contain chr21" 
    exit 1
fi

# get blacklist regions
if [[ -n ${blfilt:-} ]]; then
    bl_regions=$(cut -f 1-3 $blfilt | sed 's/\t/:/' | sed 's/\t/-/' | tr '\n' ' ')
fi

# run macs for each in array
for idx in ${!test_array[@]}; do

	test=${test_array[$idx]}
	control=${cntrl_array[$idx]:-}
	name=${name_array[$idx]}

    # check input
    if [[ $check = 'on' ]]; then
        if [[ ! -r $test ]]; then
            printf "WARNING: test BAM cannot be read: %s, skipping.\n" $test
            continue;
        elif [[ -r ${outdir}/${name}_peaks.xls ]] && [[ $force != 'on' ]]; then
            printf "WARNING: peaks already exist for %s, skipping\n" $name
            continue;
        elif [[ -n ${control:-} ]] && [[ ${control:-} != 'none' ]] && [[ ! -r ${control} ]]; then
            printf "\nERROR: control BAM cannot be read: %s, skipping.\n" $control;
            continue;
        fi
    fi

    # set basename
    test_base=$(basename "$test")
    test_base=${test_base%.*}

    # filtering options
    if [[ ${dedup} = 'on' ]]; then dedup_arg="-F 1024"; fi
    if [[ ${chrfilt} = 'on' ]]; then chr_arg=${chr_list}; fi

	# set control arg
	if [[ -n ${control:-} ]] && [[ ${control:-} != 'none' ]]; then
        ctrl_base=$(basename "$control")
        ctrl_base=${ctrl_base%.*}
        if [[ ${format} = 'BAM' ]] || [[ ${format} = 'BAMPE' ]]; then
            if [[ -n ${bl_regions:-} ]]; then
                cntrl_filt=("samtools view -b -U ${tmpdir}/${ctrl_base}.blfilt.bam"
                    "${control} ${bl_regions} > ${tmpdir}/${ctrl_base}.bl.bam;"
                    "samtools index ${tmpdir}/${ctrl_base}.blfilt.bam;"
                    "samtools view -b ${dedup_arg:-} -q ${qual}"
                    "${tmpdir}/${ctrl_base}.blfilt.bam ${chr_arg:-} >"
                    "${tmpdir}/${ctrl_base}.filt.bam;"
                    "samtools index ${tmpdir}/${ctrl_base}.filt.bam")
            else
                ctrl_filt=("samtools view -b ${dedup_arg:-} -q ${qual}"
                    "${control} ${chr_arg:-} > ${tmpdir}/${ctrl_base}.filt.bam;"
                    "samtools index ${tmpdir}/${ctrl_base}.filt.bam")
            fi
            ctrl_arg="-c ${tmpdir}/${ctrl_base}.filt.bam"
        else 
            ctrl_arg="-c ${control}"
        fi
    else
        cntrl_arg="--nolambda"
        echo "WARNING: no control - running with --nolambda."
    fi

    # set filtering commands
    if [[ ${format} = 'BAM' ]]; then
        if [[ -n ${bl_regions:-} ]]; then
            test_filt=("samtools view -b -U ${tmpdir}/${test_base}.blfilt.bam"
                "${test} ${bl_regions} > ${tmpdir}/${test_base}.bl.bam;"
                "samtools index ${tmpdir}/${test_base}.blfilt.bam;"
                "samtools view -b ${dedup_arg:-} -q ${qual}"
                "${tmpdir}/${test_base}.blfilt.bam ${chr_arg:-} >"
                "${tmpdir}/${test_base}.filt.bam;"
                "samtools index ${tmpdir}/${test_base}.filt.bam")
        else
            test_filt=("samtools view -b ${dedup_arg:-} -q ${qual}"
                "${test} ${chr_arg:-} > ${tmpdir}/${test_base}.filt.bam;"
                "samtools index ${tmpdir}/${test_base}.filt.bam")
        fi
        infile="${tmpdir}/${test_base}.filt.bam"
    else
        infile="${test}"
    fi

    # set peak calling commands
    if [[ $format != 'BEDPE' ]]; then
        predictd_cmd=("macs2 predictd"
                  "-i ${infile}"
                  "-g ${genome}"
                  "--outdir ${qcdir}"
                  "--format ${format}"
                  "--rfile ${qcdir}/${name}.predictd.R;"
                  "Rscript ${qcdir}/${name}.predictd.R")
    fi
    callpeaks_cmd=("macs2 callpeak"
                   "-t ${infile}"
                   "-g ${genome}"
                   "--outdir ${outdir}"
                   "-n ${name}"
                   "--keep-dup all"
                   "--format ${format}"
                   "${tagsize:-}"
                   "${bdg_arg:-}"
                   "${cntrl_arg:-}"
                   "${ext_arg:-}"
                   "${broad_arg:-}"
                   "${pval_arg:-}"
                   "${model_arg:-}"
                   "${shift_arg:-}"
                   "${summits_arg:-}")

    # set log file names
    logfile=${logdir}/${name}.$(basename "$0" .sh).log
    scrfile=${logdir}/${name}.$(basename "$0" .sh).scr

    # run job
    script=$(cat <<- EOS
			#!/bin/bash
			#SBATCH --time=24:00:00
			#SBATCH -N 1
			#SBATCH -n 1
			#SBATCH --mem=20G
			#SBATCH --job-name=macs2
			#SBATCH -o ${logfile}
			${depend:-}

			printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

			# load modules
			source ~/miniconda3/etc/profile.d/conda.sh
			conda activate macs2

			# check tmpdir
			mkdir -p ${tmpdir}

			# filter bams to remove dup/lowqual/scaffold
			${test_filt[@]:-}
			${ctrl_filt[@]:-}

			# run predictd for qc
			printf "\nPredicting insert size:\n"
			${predictd_cmd[@]:-}

			# call peaks
			printf "\nCalling peaks:\n"
			${callpeaks_cmd[@]}

			printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
			ls -lhAR ${tmpdir}
		EOS
	)
	echo "${script}" > ${scrfile}

	# submit job
	jobid=$(sbatch "${scrfile}")
    echo "${jobid}"
done
exit 0
