#!/bin/bash

#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err


PROGNAME=$0
usage() {
  cat << EOF >&2

Usage                               : $PROGNAME [-options] with the following options:
-h                                  : Help. Shows this text.
-d                                  : Dry run
-e <executable>                     : Path to executable (default = "")
-f <jobfile>                        : Path to simulation file, two columns formatted as [configfile seed] (default = "")
-o <output logfile>                 : Path to GNU parallel logfile (default = "")
EOF
  exit 1
}

while getopts hde:f:o: o; do
    case $o in
        (h) usage ;;
        (d) dryrun="ON";;
        (e) exec=$OPTARG;;
        (f) jobfile=$OPTARG;;
        (o) outfile=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done


if [ ! -f $jobfile ]; then
    echo "job file is not a valid file: $jobfile"
    exit 1
fi

num_cols=$(awk '{print NF}' $jobfile | head -n 1)
arg_line=$(tail -n+$SLURM_ARRAY_TASK_ID $jobfile | head -1)
config_file=$(echo $arg_line | cut -d " " -f1)
config_base=$(basename $config_file .cfg)
model_seed=$(echo $arg_line | cut -d " " -f2)
mkdir -p logs/$config_base


echo "HOSTNAME          : $HOSTNAME"
echo "CPUS ON  NODE     : $SLURM_CPUS_ON_NODE"
echo "CPUS PER NODE     : $SLURM_JOB_CPUS_PER_NODE"
echo "CPUS PER TASK     : $SLURM_CPUS_PER_TASK"
echo "MEM PER CPU       : $SLURM_MEM_PER_CPU"
echo "MEM PER NODE      : $SLURM_MEM_PER_NODE"
echo "ARRAY JOB ID      : $SLURM_ARRAY_JOB_ID"
echo "ARRAY TASK ID     : $SLURM_ARRAY_TASK_ID"
echo "JOB FILE          : $jobfile"
echo "CONFIG FILE       : $config_file"
echo "SEED              : $model_seed"


if [ "$num_cols" -eq 2 ]; then
    echo "EXEC LINE         : $exec -c $config_file -s $model_seed &> logs/$config_base/$model_seed.out"
    if [ -z  "$dryrun" ];then
      $exec -c $config_file -s $model_seed &> logs/$config_base/$model_seed.out
    fi
elif [ "$num_cols" -eq 3 ]; then
    bit_field=$(echo $arg_line | cut -d " " -f3)
    echo "EXEC LINE         : $exec -c $config_file -s $model_seed -b $bit_field &> logs/$config_base/$model_seed_$bit_field.out"
    if [ -z  "$dryrun" ];then
      $exec -c $config_file -s $model_seed -b $bit_field &> logs/$config_base/$model_seed_$bit_field.out
    fi
else
    echo "Case not implemented"
    exit 1
fi
