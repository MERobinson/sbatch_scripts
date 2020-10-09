#!/bin/bash
set -o pipefail
set -o nounset

# help message
help_message="

Simple wrapper to run BWA index

usage:
    bash $(basename "$0") [-wrfgh]
required arguments:
    -f|--fasta : path to input fasta to index
optional arguments:
    -n|--name : name prefix for index files (default = fa name)
    -o|--outdir : output directory for index files (default = pwd)
    -l|--logdir : output directory for log files (default = outdir)
    -s|--species : species name - added to sequence dictionary if provided
    -g|--genome : genome version - added to sequence dictionary if provided
additional info:
    # will generate fasta idx and dictionary if not already present

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -f|--fasta)
            fasta=$2
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
        -s|--species)
            species="-SPECIES $2"
            shift
            ;;
        -g|--genome)
            genome="-GENOME_ASSEMBLY $2"
            shift
            ;;
        *)
            printf "\nERROR: Unrecognised argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
fi

# check req files exist
if [[ ! -e ${fasta} ]]; then
	printf "\nERROR: Input fasta file does not exist: %s\n" $fasta 
	echo "$help_message"; exit 1
fi

# set/check output directories
if [[ -z ${outdir:-} ]]; then
    outdir="./"
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

# get basename/prefix
fasta_base=$(basename ${fasta})
fasta_prefix=${fasta%.fa}

# set outfile name
if [[ -z ${name:-} ]]; then
    name=${fasta_base%.fa*}
fi

# check if index & dict exist
if [[ ! -e ${fasta_prefix}.fa.fai ]]; then
    echo "WARNING: no fai file detected, will generate."
    index_cmd=("samtools faidx ${fasta}")
fi
if [[ ! -e ${fasta_prefix}.dict ]]; then
    echo "WARNING: no dict file detected, will generate."
    dict_cmd=("java -jar /opt/picard/2.21.1/picard.jar CreateSequenceDictionary"
              "REFERENCE=${fasta}"
              "OUTPUT=${fasta_prefix}.dict"
              "${species:-} ${genome:-}")
fi

# set logfiles
scr_name=$(basename $0 .sh)
scr_log="${logdir}/${name}.${scr_name}.scr"
std_log="${logdir}/${name}.${scr_name}.log"

# run job
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --walltime=24:00:00
		#SBATCH -N 1
		#SBACTH -n 1
		#SBATCH --mem=8gb
		#SBATCH --job-name=bwa_index
		#SBATCH -o ${std_log}

		# load modules
		module load SAMtools
		module load Java
		module load picard
		module load BWA
	
		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# create seq dict if req
		${dict_cmd:-}

		# create fasta index if req
		${index_cmd:-}

		# create BWA index
		bwa index ${fasta} 
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		EOS
)
echo "$script" > $scr_log

# submit job
jobid=$(qsub "$scr_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
