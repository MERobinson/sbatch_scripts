#!/bin/bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
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
if [[ ${check} = 'on' ]] && [[ ! -r ${workdir}/${transcripts} ]]; then
	printf "\nERROR: Input FASTA file does not exist: %s/%s\n" $workdir $fasta 
	echo "$help_message"; exit 1
fi

# check optional arg
if [[ ${gencode} = 'on' ]]; then
    gencode_arg='--gencode'
fi

# create required outdir
mkdir -p $workdir/$outdir

# set name
if [[ -z ${name:-} ]]; then
    name=$(basename ${transcripts})
    name="salmon_index_${name%%.*}"
fi

# set command
salmon_cmd=("salmon index -t $(basename ${transcripts})"
            "-i output/$name --type quasi ${gencode_arg:-}")

# set logfiles
scr_name=$(basename $0 .sh)
std_log=$workdir/$outdir/$name.$scr_name.std.log
pbs_log=$workdir/$outdir/$name.$scr_name.pbs.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=6gb:ncpus=1
		#PBS -j oe
		#PBS -N $name
		#PBS -q med-bio
		#PBS -o ${std_log}
	
		# load modules
		module load anaconda3/personal
		source activate salmon-0.12.0
	
		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# copy fasta to scratch
		cp -L ${workdir}/${transcripts}* .

		# run salmon
		${salmon_cmd[@]} 
		
		# copy index files to outdir
		cp -r output/${name}* $workdir/$outdir/

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
		ls -lhAR

		EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
