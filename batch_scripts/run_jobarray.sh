#!/bin/bash

#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#!#SBATCH --ntasks-per-core=1


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

mkdir -p logs
num_cols=$(awk '{print NF}' $simfile | head -n 1)
arg_line=$(tail -n+$SLURM_ARRAY_TASK_ID $simfile | head -1)


echo "Running job $SLURM_ARRAY_JOB_ID, iter $SLURM_ARRAY_TASK_ID at $HOSTNAME with inputfile $simfile"
echo "CPUS ON  NODE: $SLURM_CPUS_ON_NODE"
echo "CPUS PER NODE: $SLURM_JOB_CPUS_PER_NODE"
echo "CPUS PER TASK: $SLURM_CPUS_PER_TASK"
echo "MEM PER CPU  : $SLURM_MEM_PER_CPU"
echo "MEM PER NODE : $SLURM_MEM_PER_NODE"


if [ "$num_cols" -eq 2 ]; then
    config=$(echo $arg_line | cut -d " " -f1)
    model_seed=$(echo $arg_line | cut -d " " -f2)
    echo "Executing    : $exec -c $config -s $seed"
    $exec -c $config -s $seed
elif [ "$num_cols" -eq 3 ]; then
    config=$(echo $arg_line | cut -d " " -f1)
    model_seed=$(echo $arg_line | cut -d " " -f2)
    bit_field=$(echo $arg_line | cut -d " " -f3)
    echo "Executing    : $exec -c $config -s $seed -b $bit_field"
    $exec -c $config -s $seed -b $bit_field
else
    echo "Case not implemented"
    exit 1
fi
