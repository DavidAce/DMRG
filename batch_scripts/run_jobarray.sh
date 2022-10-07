#!/bin/bash

#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/%x-%A_%a.txt
#SBATCH --error=logs/%x-%A_%a.err

PROGNAME=$0
usage() {
  cat << EOF >&2

Usage                               : $PROGNAME [-options] with the following options:
-h                                  : Help. Shows this text.
-d                                  : Dry run
-e <executable>                     : Path to executable (default = "")
-f <jobfile>                        : Path to simulation file, two columns formatted as [configfile seed] (default = "")
-o <output logfile>                 : Path to output logfile (default = "")
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

echo "HOSTNAME                 : $HOSTNAME"
echo "JOB FILE                 : $jobfile"
echo "SLURM_CLUSTER_NAME       : $SLURM_CLUSTER_NAME"
echo "SLURM_CPUS_ON_NODE       : $SLURM_CPUS_ON_NODE"
echo "SLURM_JOB_CPUS_PER_NODE  : $SLURM_JOB_CPUS_PER_NODE"
echo "SLURM_CPUS_PER_TASK      : $SLURM_CPUS_PER_TASK"
echo "SLURM_MEM_PER_CPU        : $SLURM_MEM_PER_CPU"
echo "SLURM_MEM_PER_NODE       : $SLURM_MEM_PER_NODE"
echo "SLURM_ARRAY_JOB_ID       : $SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID      : $SLURM_ARRAY_TASK_ID"
echo "SLURM_ARRAY_TASK_STEP    : $SLURM_ARRAY_TASK_STEP"
echo "SLURM_ARRAY_TASK_MIN     : $SLURM_ARRAY_TASK_MIN"
echo "SLURM_ARRAY_TASK_MAX     : $SLURM_ARRAY_TASK_MAX"


num_cols=$(awk '{print NF}' $jobfile | head -n 1)
start_id=$SLURM_ARRAY_TASK_ID
end_id=$((SLURM_ARRAY_TASK_ID + SLURM_ARRAY_TASK_STEP - 1))
exit_code_save=0

echo "TASK ID SEQUENCE         : $(seq -s ' ' $start_id $end_id)"

for id in $(seq $start_id $end_id); do
  arg_line=$(tail -n+$id $jobfile | head -1)
  if [ -z "$arg_line" ]; then
    continue
  fi
  config_file=$(echo "$arg_line" | cut -d " " -f1)
  config_base=$(echo "$config_file" | xargs -l basename)
  config_dir=$(echo "$config_file" | xargs -l dirname)
  model_seed=$(echo "$arg_line" | cut -d " " -f2)
  logdir=logs/$config_dir/$config_base
  mkdir -p $logdir
  echo "JOB ID                   : $id"
  echo "TIME                     : $(/bin/date)"
  echo "CONFIG LINE              : $arg_line"
  if [ "$num_cols" -eq 2 ]; then
      echo "EXEC LINE                : $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed &>> $logdir/$model_seed.txt"
      if [ -z  "$dryrun" ];then
        $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed &>> $logdir/$model_seed.txt
        exit_code=$?
        echo "EXIT CODE                : $exit_code"
        if [ "$exit_code" != "0" ]; then
          exit_code_save=$exit_code
          continue
        fi
      fi
  elif [ "$num_cols" -eq 3 ]; then
      bit_field=$(echo $arg_line | cut -d " " -f3)
      echo "BITFIELD                 : $bit_field"
      echo "EXEC LINE                : $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed -b $bit_field &>> $logdir/$model_seed_$bit_field.txt"
      if [ -z  "$dryrun" ];then
        $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed -b $bit_field &>> $logdir/$model_seed_$bit_field.txt
        exit_code=$?
        echo "EXIT CODE         : $exit_code"
        if [ "$exit_code" != "0" ]; then
          exit_code_save=$exit_code
          continue
        fi
      fi
  else
      echo "Case not implemented"
      exit 1
  fi
done

exit $exit_code_save


