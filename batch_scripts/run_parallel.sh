#!/bin/bash

#SBATCH --nodes=1
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/%x-%A.out
#SBATCH --error=logs/%x-%A.err




PROGNAME=$0
usage() {
  cat << EOF >&2

Usage                               : $PROGNAME [-options] with the following options:
-h                                  : Help. Shows this text.
-e <executable>                     : Path to executable (default = "")
-f <simfile>                        : Path to simulation file, two columns formatted as [configfile seed] (default = "")
-o <output logfile>                 : Path to GNU parallel logfile (default = "")
EOF
  exit 1
}

while getopts he:f:o: o; do
    case $o in
        (h) usage ;;
        (e) exec=$OPTARG;;
        (f) simfile=$OPTARG;;
        (o) outfile=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done


if [ ! -f $simfile ]; then
    echo "Simfile is not a valid file: $simfile"
    exit 1
fi

if [ -z "$outfile" ] ; then
    simbase=$(basename $simfile .sim)
    outfile=logs/$simbase.log
fi



echo "Running job $SLURM_JOB_ID at $HOSTNAME with simfile $simfile"
outdir=$(dirname $outfile)
mkdir -p $outdir/$simbase


num_cols=$(awk '{print NF}' $simfile | head -n 1)

if [ "$num_cols" -eq 2 ]; then
    cat $simfile | parallel --memfree $SLURM_MEM_PER_CPU --joblog $outfile --colsep ' ' "$exec -i {1} -r {2} &> $outdir/${simbase}/{1./}_{2}.out"
elif [ "$num_cols" -eq 3 ]; then
    cat $simfile | parallel --memfree $SLURM_MEM_PER_CPU --joblog $outfile --colsep ' ' "$exec -i {1} -r {2} -s {3} &> $outdir/${simbase}/{1./}_{2}_{3}.out"
else
    echo "Case not implemented"
    exit 1
fi