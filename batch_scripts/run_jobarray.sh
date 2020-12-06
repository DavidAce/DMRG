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

echo "HOSTNAME          : $HOSTNAME"
echo "CLUSTER           : $SLURM_CLUSTER_NAME"
echo "CPUS ON  NODE     : $SLURM_CPUS_ON_NODE"
echo "CPUS PER NODE     : $SLURM_JOB_CPUS_PER_NODE"
echo "CPUS PER TASK     : $SLURM_CPUS_PER_TASK"
echo "MEM PER CPU       : $SLURM_MEM_PER_CPU"
echo "MEM PER NODE      : $SLURM_MEM_PER_NODE"
echo "ARRAY JOB ID      : $SLURM_ARRAY_JOB_ID"
echo "ARRAY TASK ID     : $SLURM_ARRAY_TASK_ID"
echo "ARRAY TASK STEP   : $SLURM_ARRAY_TASK_STEP"
echo "ARRAY TASK MIN ID : $SLURM_ARRAY_TASK_MIN"
echo "ARRAY TASK MAX ID : $SLURM_ARRAY_TASK_MAX"
echo "JOB FILE          : $jobfile"


num_cols=$(awk '{print NF}' $jobfile | head -n 1)
start_id=$SLURM_ARRAY_TASK_ID
end_id=$((SLURM_ARRAY_TASK_ID+SLURM_ARRAY_TASK_STEP))
exit_code_save=0
for id in $(seq $start_id $end_id); do
  arg_line=$(tail -n+$id $jobfile | head -1)
  config_file=$(echo "$arg_line" | cut -d " " -f1)
  config_base=$(echo "$config_file" | xargs -l basename)
  config_dir=$(echo "$config_file" | xargs -l dirname)
  model_seed=$(echo "$arg_line" | cut -d " " -f2)
  logdir=logs/$config_dir/$config_base
  mkdir -p $logdir

  echo "JOB FILE LINE(S)  : $arg_line"
  echo "CONFIG FILE       : $config_file"
  echo "SEED              : $model_seed"
  if [ "$num_cols" -eq 2 ]; then
      echo "EXEC LINE         : $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed &>> $logdir/$model_seed.out"
      if [ -z  "$dryrun" ];then
        $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed &>> $logdir/$model_seed.out
        exit_code=$?
        if [ "$exit_code" == "0" ]; then
          echo "SUCCESS           : $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed &>> $logdir/$model_seed.out"
        else
          echo "EXIT CODE         : $exit_code"
          echo "FAILED            : $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed &>> $logdir/$model_seed.out"
          exit_code_save=$exit_code
          continue
        fi
      fi
  elif [ "$num_cols" -eq 3 ]; then
      bit_field=$(echo $arg_line | cut -d " " -f3)
      echo "EXEC LINE         : $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed -b $bit_field &>> $logdir/$model_seed_$bit_field.out"
      if [ -z  "$dryrun" ];then
        $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed -b $bit_field &>> $logdir/$model_seed_$bit_field.out
        exit_code=$?
        if [ "$exit_code" == "0" ]; then
          echo "SUCCESS           : $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed -b $bit_field &>> $logdir/$model_seed_$bit_field.out"
        else
          echo "EXIT CODE         : $exit_code"
          echo "FAILED            : $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed -b $bit_field &>> $logdir/$model_seed_$bit_field.out"
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


