#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
check='on'
short='off'
scr_name=$(basename "$0")

# help message
help_message="
Wrapper to align reads with BWA

usage:
    bash $scr_name [-options] -fq1 <FASTQ>
required arguments:
    -fq1|--fastq1 : FASTQ filename of read 1 
    -F|--fasta : reference FASTA used in index generation
    -I|--index : BWA index prefix to align against 
optional arguments:
    -fq2|--fastq2 : FASTQ filename of read 2 (if paired end)
    -s|--short : align in short read SE mode [on|off] (default = off)
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on,off] (default = on)
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
        -s|--short)
            short=$2
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
        printf "\nERROR: FASTQ r1 cannot be read: %s\n" $fq1
        echo "$help_message"; exit 1
    elif [[ -n "${fq2:-}" ]] && [[ ! -r "${fq2:-}" ]]; then
        printf "\nERROR: FASTQ r2 cannot be read: %s\n" ${fq2:-}
        echo "$help_message"; exit 1
    elif [[ ! -r $fasta ]]; then
        printf "\nERROR: FASTA file cannot be read: %s\n" $fasta
        echo "$help_message"; exit 1
    elif [[ ! -e $index.bwt ]]; then
        printf "\nERROR: Index file does not exist: %s\n" $index
        echo "$help_message"; exit 1
    elif [[ ! -e "${fasta%.fa*}.dict" ]]; then
        printf "\nERROR: No FASTA dict found next to FASTA: %s\n" $fasta
    fi
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

# set basenames
fq1_base=$(basename "$fq1")

# set PE/SE arguments
if [[ -n "${fq2:-}" ]]; then
    fq2_fq2sam="FASTQ2=${fq2}"
    fq2_bwamem="-p"
fi

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    echo "WARNING: name not set, extracting from: $fq1"
    name=${fq1_base%%.*}
fi

# set sample name to name if not provided
if [[ -z "${sm:-}" ]]; then
    sm="SAMPLE_NAME=$name"
fi

# extract read group info (if not provided)
if [[ ${pl:-} = "PLATFORM=ILLUMINA" ]]; then
    ext=${fq1##*.}
    if [[ ${ext} == "gz"  ]]; then
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

# set commands
fastqc_command=("fastqc --noextract --dir ${tmpdir} -o ${qcdir} -t 20 ${fq1} ${fq2:-}")
fq2sam_command=("java -Xmx32G -jar /opt/picard/2.21.1/picard.jar FastqToSam" 
                "FASTQ=${fq1}" 
                "OUTPUT=${tmpdir}/${name}.bam"
                "TMP_DIR=${tmpdir}"
                "${fq2_fq2sam:-} ${id:-} ${lb:-} ${pl:-}" 
                "${pu:-} ${pm:-} ${cn:-} ${dt:-} ${pi:-} ${sm:-}")
madapt_command=("java -Xmx32G -jar /opt/picard/2.21.1/picard.jar MarkIlluminaAdapters"
                "I=${tmpdir}/${name}.bam"
                "O=${tmpdir}/${name}.madapt.bam"
                "M=${qcdir}/${name}.mark_adapters_metrics"
                "TMP_DIR=${tmpdir}")
sam2fq_command=("java -Xmx32G -jar /opt/picard/2.21.1/picard.jar SamToFastq"
                "I=${tmpdir}/${name}.madapt.bam" 
                "FASTQ=${tmpdir}/${name}.unaligned.fq" 
                "CLIPPING_ATTRIBUTE=XT" 
                "CLIPPING_ACTION=2" 
                "INTERLEAVE=true" 
                "NON_PF=true"
                "TMP_DIR=${tmpdir}")
if [[ $short = 'on' ]]; then
    sai_command=("bwa aln ${index} ${tmpdir}/${name}.unaligned.fq"
                 "> ${tmpdir}/${name}.unaligned.sai")
    align_command=("bwa samse ${index} ${tmpdir}/${name}.unaligned.sai"
                   "${tmpdir}/${name}.unaligned.fq > ${tmpdir}/${name}.aligned.sam")
else
    align_command=("bwa mem -M -t 20 ${fq2_bwamem:-} ${index}"
                   "${tmpdir}/${name}.unaligned.fq > ${tmpdir}/${name}.aligned.sam")
fi
merge_command=("java -Xmx32G -jar /opt/picard/2.21.1/picard.jar MergeBamAlignment"
                "R=${fasta}" 
                "UNMAPPED_BAM=${tmpdir}/${name}.madapt.bam" 
                "ALIGNED=${tmpdir}/${name}.aligned.sam"
                "O=${tmpdir}/${name}.bam"
                "SORT_ORDER=coordinate"
                "CREATE_INDEX=true"
                "ADD_MATE_CIGAR=true"
                "CLIP_ADAPTERS=false"
                "CLIP_OVERLAPPING_READS=true"
                "INCLUDE_SECONDARY_ALIGNMENTS=true" 
                "MAX_INSERTIONS_OR_DELETIONS=-1"
                "PRIMARY_ALIGNMENT_STRATEGY=MostDistant" 
                "ATTRIBUTES_TO_RETAIN=XS"
                "TMP_DIR=${tmpdir}")
mdup_command=("sambamba markdup" 
            	"--nthreads=20"
            	"--tmpdir=${tmpdir}"
                "${tmpdir}/${name}.bam ${outdir}/${name}.bam")
alstat_command=("java -Xmx16g -jar /opt/picard/2.21.1/picard.jar CollectAlignmentSummaryMetrics"
            	"R=${fasta}"
            	"I=${outdir}/${name}.bam"
            	"O=${qcdir}/${name}.alignment_summary_metrics"
                "TMP_DIR=${tmpdir}")


# set log file names
logfile=${logdir}/${name}.$(basename "$0" .sh).log
scrfile=${logdir}/${name}.$(basename "$0" .sh).scr

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH -n 20
		#SBATCH --mem=32gb
		#SBATCH --job-name=$name.bwa
		#SBATCH -o ${logfile}
		${depend:-}

		# load modules
		#module load FastQC
		#module load SAMtools
		#module load Java
		#module load picard
		#module load BWA
		#module load Sambamba
		source ~/miniconda3/etc/profile.d/conda.sh
		conda activate bwa

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		mkdir -p ${tmpdir}

		# run fastqc
		printf "\nRunning FASTQC:\n" \`date '+%H:%M:%S'\` 
		${fastqc_command[@]}

		# covert fastq to ubam
		printf "\nConverting to uBAM:\n" \`date '+%H:%M:%S'\` 
		${fq2sam_command[@]}
		
		# mark adapters
		printf "\nMarking adapters:\n" \`date '+%H:%M:%S'\`
		${madapt_command[@]}

		# convert uBAM to interleaved fastq
		printf "\nConverting to FASTQ:\n" \`date '+%H:%M:%S'\`
		${sam2fq_command[@]}

		# align to genome with bwa
		printf "\nAligning to genome:\n" \`date '+%H:%M:%S'\`
		${sai_command[@]:-}
		${align_command[@]}

		# merge uBAM and aligned
		printf "\nMerging aligned with uBAM:\n" \`date '+%H:%M:%S'\`
		${merge_command[@]}

		# # mark dup
		# printf "\nMarking duplicates:\n" \`date '+%H:%M:%S'\`
		# ${mdup_command[@]}

		# sort & index final bam & copy to outdir
		samtools sort -@ 12 -o ${outdir}/${name}.bam ${tmpdir}/${name}.bam
		samtools index ${outdir}/${name}.bam 

		# alignment metrics
		printf "\nCollecting alignment metrics:\n" \`date '+%H:%M:%S'\`
		${alstat_command[@]} 
 
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
        
	EOS
) 
echo "$script" > $scrfile

# submit job
jobid=$(sbatch "$scrfile")

# echo job id and exit
echo "JOBID: $jobid"
exit 0 
