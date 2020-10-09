#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
check='on'
clip_act='X'
scr_name=$(basename "$0")

# help message
help_message="
Wrapper to align reads with Bowtie2

usage:
    bash $scr_name [-options] -fq1 <FASTQ> -I <INDEX> -F <FASTA>
required arguments:
    -fq1|--fastq1 : FASTQ filename of read 1 
    -F|--fasta : reference FASTA used in index generation
    -I|--index : BWA index prefix to align against 
optional arguments:
    -fq2|--fastq2 : FASTQ filename of read 2 (if paired end)
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --bt2_arg : options to pass to bowtie2
    --clip_act : clipping action setting for Picard SamToFastq [N|X|INT] (default = X)
    --check : whether to check input files [on|off] (default = on)
    --depend : list of PBS dependencies (default = NULL)
bam header arguments:
    --sm : sample name (default = name)
    --id : read group ID (usually flow cell + lane)
    --lb : library name (name unique to each library)
    --pu : platform unit (usually flow cell + barcode + lane)
    --pl : sequencing platform sequenced on, valid options:
           [ILLUMINA,PACBIO,IONTORRENT,ONT,CAPILLARY,LS454,SOLID,HELICOS]
    --cn : sequencing centre sequencing was performed at 
    --dt : run date (Iso8601Date)
    --pi : predicted median insert size (e.g. 200)
    --pm : platform model - further discription of platform
additional info:
    # bt2_arg should be a quoted string of bowtie2 options, e.g. '-k 4 -X 500 -I 50' 
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified
    # any of the bam header arguments provided will be included in the header
    # if platform is ILLUMINA - will automatically extract PU and ID from read name 

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -fq1|--fastq1)
            fq1=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -I|--index)
            index=$2
            shift
            ;;
        -fq2|--fastq2)
            fq2=$2
            shift
            ;;
        -o|--outdir)
            outdir="./$2"
            shift
            ;;
        -q|--qcdir)
            qcdir="./$2"
            shift
            ;;
        -l|--logdir)
            logdir="./$2"
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        --clip_act)
            clip_act=$2
            shift
            ;;
        --bt2_arg)
            bt2_arg=$2
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
        --sm)
            sm="SAMPLE_NAME=$2"
            shift
            ;;
        --id)
            id="READ_GROUP_NAME=$2"
            shift
            ;;
        --cn)
            cn="SEQUENCING_CENTER=$2"
            shift
            ;;
        --dt)
            dt="RUN_DATE=$2"
            shift
            ;;
        --lb)
            lb="LIBRARY_NAME=$2"
            shift
            ;;
        --pi)
            pi="PREDICTED_INSERT_SIZE=$2"
            shift
            ;;
        --pl)
            pl="PLATFORM=$2"
            shift
            ;;
        --pm)
            pm="PLATFORM_MODEL=$2"
            shift
            ;;
        --pu)
            pu="PLATFORM_UNIT=$2"
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
if [[ -z ${fq1:-} ]]; then
    printf "\nERROR: --fastq1 argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${index:-} ]]; then
    printf "\nERROR: --index argument required\n"
    echo "$help_message"; exit 1
fi

# check files 
if [[ "${check:-}" = on ]]; then
    if [[ ! -r $fq1 ]]; then
        printf "\nERROR: FASTQ r1 cannot be read:%s\n" $fq1
        echo "$help_message"; exit 1
    elif [[ -n "${fq2:-}" ]] && [[ ! -r "${fq2:-}" ]]; then
        printf "\nERROR: FASTQ r2 cannot be read: %s\n" ${fq2:-}
        echo "$help_message"; exit 1
    elif [[ ! -r $fasta ]]; then
        printf "\nERROR: FASTA file cannot be read: %s\n" $fasta
        echo "$help_message"; exit 1
    elif [[ ! -e $index.1.bt2 ]]; then
        printf "\nERROR: Index file does not exist: %s\n" $index
        echo "$help_message"; exit 1
    fi
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
if [[ -z "${qcdir:-}" ]]; then
    qcdir=${outdir}
fi
mkdir -p $qcdir
tmpdir=$(mktemp -d)

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    name=$(basename $fq1)
    name=${name%%.*}
    echo "WARNING: name not set, extracted from FASTQ: $name"
fi

# set optional arguments
if [[ -n "${fq2:-}" ]]; then
    fq2_trim_out="-p ${tmpdir}/${name}.2.fastq"
    fq2_adapt="-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    fq2_fq2sam="FASTQ2=${tmpdir}/${name}.2.fastq"
    bt2_input="--interleaved ${tmpdir}/${name}.unaligned.fastq"
else
    bt2_input="-U ${tmpdir}/${name}.unaligned.fastq"
fi

# set sample name to name if not provided
if [[ -z "${sm:-}" ]]; then
    sm="SAMPLE_NAME=$name"
fi

# extract read group info (if not provided)
if [[ ${pl:-} = "PLATFORM=ILLUMINA" ]]; then
    ext=${fq1##*.}
    if [[ $ext == "gz" ]]; then 
        IFS=":" read -r -a read_name <<< $(gzip -dc $fq1 | head -n 1)
    else
        IFS=":" read -r -a read_name <<< $(head -n 1 $fq1)
    fi
    fcid=${read_name[2]} # flow cell ID
    fcln=${read_name[3]} # flow cell lane
    ridx=${read_name[9]} # barcode/index
    if [[ -z ${id:-} ]]; then id="READ_GROUP_NAME=$fcid.$fcln"; fi
    if [[ -z ${pu:-} ]]; then pu="PLATFORM_UNIT=$fcid.$fcln.$ridx"; fi
fi

# deal with memory issue - temp
export _JAVA_OPTIONS=-Djava.io.tmpdir=${tmpdir}

# set commands
picard_path="~/miniconda3/envs/bowtie2/share/picard-2.21.8-0"
fastqc_command=("fastqc --noextract --dir ${tmpdir} -o ${qcdir} -t 18 ${fq1} ${fq2:-}")
cutadapt_cmd=("cutadapt --minimum-length 15 -e 0.1 --trim-n -j 12"
              "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA ${fq2_adapt:-}"
              "-o ${tmpdir}/${name}.1.fastq ${fq2_trim_out:-} ${fq1} ${fq2:-}")
fq2sam_command=("java -Xmx32G -jar ${picard_path}/picard.jar FastqToSam" 
                "FASTQ=${tmpdir}/${name}.1.fastq" 
                "OUTPUT=${tmpdir}/${name}.bam"
                "TMP_DIR=${tmpdir}"
                "${fq2_fq2sam:-} ${id:-} ${lb:-} ${pl:-}" 
                "${pu:-} ${pm:-} ${cn:-} ${dt:-} ${pi:-} ${sm:-}")
madapt_command=("java -Xmx32G -jar ${picard_path}/picard.jar MarkIlluminaAdapters"
                "I=${tmpdir}/${name}.bam"
                "O=${tmpdir}/${name}.madapt.bam"
                "M=${qcdir}/${name}.mark_adapters_metrics")
sam2fq_command=("java -Xmx32G -jar ${picard_path}/picard.jar SamToFastq"
                "I=${tmpdir}/${name}.madapt.bam" 
                "FASTQ=${tmpdir}/${name}.unaligned.fastq" 
                "CLIPPING_ATTRIBUTE=XT" 
                "CLIPPING_ACTION=${clip_act}" 
                "CLIPPING_MIN_LENGTH=20"
                "INTERLEAVE=true"
                "NON_PF=true"
                "TMP_DIR=${tmpdir}")
align_command=("bowtie2 ${bt2_arg:-}"
                "--threads 18 -q" 
                "-x ${index}" 
                "${bt2_input}"
                "-S ${tmpdir}/${name}.aligned.sam")
merge_command=("java -Xmx32G -jar ${picard_path}/picard.jar MergeBamAlignment"
                "R=${fasta}" 
                "UNMAPPED_BAM=${tmpdir}/${name}.madapt.bam" 
                "ALIGNED=${tmpdir}/${name}.aligned.sort.bam"
                "O=${tmpdir}/${name}.bam"
                "SORT_ORDER=coordinate"
                "CREATE_INDEX=true"
                "ADD_MATE_CIGAR=true"
                "CLIP_ADAPTERS=true"
                "CLIP_OVERLAPPING_READS=true"
                "INCLUDE_SECONDARY_ALIGNMENTS=true" 
                "MAX_INSERTIONS_OR_DELETIONS=-1"
                "PRIMARY_ALIGNMENT_STRATEGY=MostDistant" 
                "ATTRIBUTES_TO_RETAIN=XS"
                "MIN_UNCLIPPED_BASES=20"
                "UNMAP_CONTAMINANT_READS=true"
                "TMP_DIR=${tmpdir}")
mdup_command=("sambamba markdup" 
            	"--nthreads=18"
            	"--tmpdir=${tmpdir}"
                "${tmpdir}/${name}.bam"
                "${tmpdir}/${name}.mdup.bam")
alstat_command=("java -Xmx16g -jar ${picard_path}/picard.jar CollectAlignmentSummaryMetrics"
            	"R=${fasta}"
            	"I=${outdir}/${name}.bam"
            	"O=${qcdir}/${name}.alignment_summary_metrics"
                "TMP_DIR=${tmpdir}")

# set log file names
scr_name=$(basename ${0} .sh)
std_log="${logdir}/${name}.${scr_name}.log"
scr_log="${logdir}/${name}.${scr_name}.scr"

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 18
		#SBATCH --mem=24gb
		#SBATCH --job-name=${name}.bt2
		#SBATCH -o ${std_log}
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		mkdir -p ${tmpdir}

		# load env
		source ~/miniconda3/etc/profile.d/conda.sh
		conda activate bowtie2

		# run fastqc
		printf "\nRunning FASTQC:\n"  
		${fastqc_command[@]}

		# trim adapters
		printf "\nTrimming adapters:\n"
		${cutadapt_cmd[@]} > ${qcdir}/${name}.cutadapt_metrics.txt

		# covert fastq to ubam
		printf "\nConverting to uBAM:\n"
		${fq2sam_command[@]}

		# mark adapters
		printf "\nMarking adapters:\n"
		${madapt_command[@]}

		# convert uBAM to interleaved fastq
		printf "\nConverting to FASTQ:\n"
		${sam2fq_command[@]}

		# align to genome with bwa
		printf "\nAligning to genome:\n"
		${align_command[@]}

		samtools sort -n -O bam -o ${tmpdir}/${name}.aligned.sort.bam \
			${tmpdir}/${name}.aligned.sam

		# merge uBAM and aligned
		printf "\nMerging aligned with uBAM:\n"
		${merge_command[@]}

		# mark dup
		printf "\nMarking duplicates:\n"
		${mdup_command[@]}

		# sort & index final bam & copy to outdir
		samtools sort -o ${outdir}/${name}.bam ${tmpdir}/${name}.mdup.bam
		samtools index ${outdir}/${name}.bam 

		# alignment metrics
		printf "\nCollecting alignment metrics:\n"
		${alstat_command[@]} 
 
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
		ls -lhAR ${tmpdir}
        rm -r ${tmpdir}
	EOS
) 
echo "$script" > $scr_log

# submit job
jobid=$(sbatch "$scr_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0 
