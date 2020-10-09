#!/bin/bash
set -o pipefail
set -o nounset

# default arg
outdir=$PWD
gencode='on'
check='on'

# help message
help_message="

wrapper to generate Salmon-0.12.0 index.

usage:
    bash $(basename "$0") [-options] -t <transcripts>
required arguments:
    -t|--transcripts : transcripts [FASTA]
optional arguments:
    -n|--name : name prefix for index files (default = fa name)
    -o|--outdir : output directory for index files (default = pwd)
    --gencode : whether input fasta is gencode format [on|off] (default = on)
    --check : whether to check input files [on|off] (default = on)

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -t|--transcripts)
            transcripts=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        --gencode)
            gencode=$2
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        *)
            echo "Error: Illegal argument"
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z ${transcripts:-} ]]; then
    printf "\nERROR: --transcripts argument required\n"
    echo "$help_message"; exit 1
fi

# check input files exist
if [[ ${check} = 'on' ]] && [[ ! -r ${transcripts} ]]; then
	printf "\nERROR: Input FASTA file does not exist: %s\n" $fasta 
	echo "$help_message"; exit 1
fi

# check optional arg
if [[ ${gencode} = 'on' ]]; then
    gencode_arg='--gencode'
fi

# create outdir
mkdir -p $outdir

# set name
if [[ -z ${name:-} ]]; then
    name=$(basename ${transcripts})
    name="salmon_index_${name%%.*}"
fi

# set command
salmon_cmd=("salmon index -p 12 ${gencode_arg:-}"
            "-t ${transcripts} -i ${outdir}/$name")

# set logfiles
scr_name=$(basename $0 .sh)
std_log=$outdir/$name.$scr_name.std.log
scr_log=$outdir/$name.$scr_name.scr.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#SBATCH --time=24:00:00
		#SBATCH -N 1
		#SBATCH --mem=6gb
		#SBATCH -n 12
		#SBATCH --job-name=salmon_index
	
		# load modules
		source ~/miniconda3/etc/profile.d/conda.sh
		conda activate salmon
	
		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# run salmon
		${salmon_cmd[@]} 
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		EOS
)
echo "$script" > $scr_log

# submit job
jobid=$(qsub "$scr_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
