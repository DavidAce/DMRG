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
config_file=$(echo $arg_line | cut -d " " -f1)
config_base=$(basename $config_file .cfg)
model_seed=$(echo $arg_line | cut -d " " -f2)

echo "HOSTNAME          : $HOSTNAME"
echo "CPUS ON  NODE     : $SLURM_CPUS_ON_NODE"
echo "CPUS PER NODE     : $SLURM_JOB_CPUS_PER_NODE"
echo "CPUS PER TASK     : $SLURM_CPUS_PER_TASK"
echo "MEM PER CPU       : $SLURM_MEM_PER_CPU"
echo "MEM PER NODE      : $SLURM_MEM_PER_NODE"
echo "ARRAY JOB ID      : $SLURM_ARRAY_JOB_ID"
echo "ARRAY TASK ID     : $SLURM_ARRAY_TASK_ID"
echo "INPUT FILE        : $simfile"
echo "CONFIG FILE       : $config_file"
echo "SEED              : $model_seed"

if [ "$num_cols" -eq 2 ]; then
    echo "EXEC LINE         : $exec -c $config_file -s $model_seed &> logs/$config_base/$model_seed.out"
    $exec -c $config_file -s $model_seed &> logs/$config_base/$model_seed.out
elif [ "$num_cols" -eq 3 ]; then
    echo "EXEC LINE         : $exec -c $config_file -s $model_seed -b $bit_field &> logs/$config_base/$model_seed_$bit_field.out"
    bit_field=$(echo $arg_line | cut -d " " -f3)
    $exec -c $config_file -s $model_seed -b $bit_field &> logs/$config_base/$model_seed_$bit_field.out
else
    echo "Case not implemented"
    exit 1
fi
