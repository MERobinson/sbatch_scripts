#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
outdir=''
bamdir='bam'
logdir='logs'
qcdir='qc/metrics'
fqdir='fastq'
scriptdir='scripts'
trackdir='tracks'
sradir='sradata'

# help message
help_message="

Pipeline to fetch SRA, extract FASTQ, align with BWA MEM and generate tracks.

usage:
    bash $(basename "$0") [-options] -s <SRA_accession>
required arguments:
    -s|--sra : comma separated list of SRA run accession numbers
    -F|--fasta : whole genome fasta file 
    -I|--index : BWA index file base for alignment
    -g|--genome : genome version [hg19,hg38,mm9,mm10]
optional arguments:
    -bd|--bamdir : output directory for bam files (default = bam)
    -ld|--logdir : output directory for log files (default = logs)
    -fd|--fqdir : output directory for fastq files (default = fastq)
    -qd|--qcdir : output directory for qc files (default = qc)
    -sd|--scriptdir : directory containing scripts required (default = scripts)
    -td|--trackdir : output directory for genome tracks (default = tracks)
    -n|--name : comma separated list of names for output files (default = SRA accession)
    -e|--extend : bp extension for generating coverage track (default = 200)
additional info:
    # all paths should be relative to working directory
    # required scripts = fetch_fastq_sra.sh, align_fastq_bwa.sh & generate_tracks_deeptools.sh 
    # either SRA run IDs or GSE sample IDs of the form SRRxxxxxx or GSMxxxxxx can be supplied
    # run info will be retrieved for all runs of a given sample and merged post-alignment

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -s|--sra)
            sra_list=$2
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
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -bd|--bamdir)
            bamdir=$2
            shift
            ;;
        -ld|--logdir)
            logdir=$2
            shift
            ;;
        -qd|--qcdir)
            qcdir=$2
            shift
            ;;
        -fd|--fqdir)
            fqdir=$2
            shift
            ;;
        -sd|--scriptdir)
            scriptdir=$2
            shift
            ;;
        -td|--trackdir)
            trackdir=$2
            shift
            ;;
        -g|--genome)
            genome=$2
            shift
            ;;
        -n|--name)
            names_list=$2
            shift
            ;;
        -e|--extend)
            extend=$2
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check/parse required arguments
if [[ -z ${sra_list:-} ]]; then
    printf "\nERROR: --sra argument required.\n"
    echo "$help_message"; exit 1
else
    IFS=',' read -r -a sra_array <<< "$sra_list"
fi
if [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required.\n"
    echo "$help_message"; exit 1
elif [[ ! -r $fasta ]]; then
    printf "\nERROR: Input FASTA file cannot be read: %s\n" $fasta
    echo "$help_message"; exit 1
fi
if [[ -z ${index:-} ]]; then
    printf "\nERROR: --index argument required.\n"
    echo "$help_message"; exit 1
elif [[ ! -r "${index}.ann" ]]; then
    printf "\nERROR: Input index files cannot be read: %s\n" $index
    echo "$help_message"; exit 1
fi

# check sample names
if [[ -z ${names_list:-} ]]; then
	names_array=(${sra_array[@]})
    echo "WARNING: No name argument given, using accession numbers."
else
	IFS=',' read -r -a names_array <<< "$names_list"
fi

# check genome
if [[ -n ${genome:-} ]]; then
    genome_arg="--genome $genome"
fi

# create required dirs
mkdir -p $logdir
mkdir -p $qcdir
mkdir -p $bamdir
mkdir -p $trackdir

# load sra toolkit
module load SRA-Toolkit
module load edirect

# for each SRA number, submit PBS job
for idx in ${!sra_array[@]}; do

    unset merge_depend
    unset merge_arg
    unset depend

	sra=${sra_array[$idx]}
	name=${names_array[$idx]}
    printf "\nProcessing sample %s\n" $name
 
    # add genome version to name if given
    if [[ -n $genome ]]; then
        comb_name=$name.$genome
    else
        comb_name=$name
    fi

	# fetch SRA run info
    if [[ ! -r $logdir/$name.run_info.csv ]]; then
        $(esearch -db sra -query $sra | efetch -format runinfo \
          > $logdir/$name.run_info.csv)
    fi

    # check n runs
    n_runs=$(wc -l $logdir/$name.run_info.csv | cut -f1 -d' ') 
    if [[ $n_runs < 2 ]]; then
        echo "WARNING: no runs found for sample $sra"
        continue
    fi

    if [[ ! -r $bamdir/$comb_name.bam ]]; then

        # fetch and align per run
        while read line; do

            unset depend

            IFS="," read -r -a srr_info <<< "$line"
            echo $line
            srr=$(echo "$line" | grep -Eo "^SRR[0-9]+")
            libtype=$(echo "$line" | grep -Eo "(SINGLE|PAIRED)")
            platform=$(echo "$line" | grep -Eo "ILLUMINA")
            run_name="${name}_${srr}"

            if [[ $libtype = PAIRED ]]; then
                fq2_arg="-fq2 $fqdir/${run_name}_2.fastq.gz"
            fi

            printf "\tprocessing run: %s\n" $srr
            printf "\tLIBTYPE: %s\n" $libtype

            # get fastq
            if [[ -r $fqdir/${run_name}_1.fastq.gz ]]; then
                printf "\tFASTQ already exists, skipping fetch\n"
            else
                sra_call=("bash $scriptdir/fetch_fastq_sra.sh"
                          "--sra $srr --outdir $fqdir"
                          "--logdir $logdir --name $run_name")
                ${sra_call[@]} &> tmp.log
                if [[ $? -ne 0 ]]; then
                    printf "\nERROR: FASTQ fetching failed for sample: %s\n" $srr
                    echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
                else
                    jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
                    depend="afterok:$jobid"
                fi
            fi

            # align
            if [[ -r $bamdir/${run_name}.bam ]]; then
                printf "\tRun BAM file already exists, skipping alignment\n"
            else
                bwa_call=("bash ${scriptdir}/align_fastq_bwa.sh"
                          "-fq1 ${fqdir}/${run_name}_1.fastq.gz ${fq2_arg:-}"
                          "--fasta ${fasta} --index ${index}"
                          "--outdir ${bamdir} --qcdir ${qcdir}"
                          "--logdir ${logdir} --check no --name ${run_name}")
                if [[ -n ${depend:-} ]]; then bwa_call="${bwa_call} --depend $depend"; fi
                ${bwa_call[@]} &> tmp.log
                if [[ $? -ne 0 ]]; then
                    printf "\nERROR: Alignment failed for sample: %s\n" $run_name
                    echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
                else
                    jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
                    merge_depend="${merge_depend:-},afterok:$jobid"
                fi
            fi
            merge_arg="${merge_arg:-},$bamdir/${run_name}.bam"

        done <<< "$(tail -n +2 $logdir/$name.run_info.csv)"

        # merge runs
        merge_arg=${merge_arg#*,}
        if [[ -n ${merge_depend:-} ]]; then
            merge_depend="--depend ${merge_depend#*,} --check no"
        fi 
        merge_call=("bash $scriptdir/merge_bam_picard.sh"
                    "--bam_list $merge_arg --fasta $fasta"
                    "--name $comb_name --outdir $bamdir"
                    "--qcdir $qcdir --logdir $logdir ${merge_depend:-}")
        ${merge_call[@]} &> tmp.log
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: Merging failed for sample: $name\n"
            echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            depend="afterok:$jobid"
        fi
    else
        printf "\tFinal BAM already exists for $sra\n"
    fi 

    # generate track
    if [[ -r $trackdir/$comb_name.bedGraph ]]; then
        printf "\tTrack already exists, skipping pileup\n"
    else
        if [[ -n ${depend:-} ]]; then
            depend_arg="--depend $depend --check no"
        fi
        track_call=("bash $scriptdir/generate_track_deeptools.sh"
                    "--bam $bamdir/$comb_name.bam --outdir $trackdir"
                    "--fasta $fasta --logdir $logdir ${depend_arg:-}"
                    "--name $comb_name ${genome_arg:-}")
        ${track_call[@]} &> tmp.log
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: Track generation failed for sample: %s\n" $sra
            echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
        fi
    fi

    if [[ -r tmp.log ]]; then rm tmp.log; fi

done

