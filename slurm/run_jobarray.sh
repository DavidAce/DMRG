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
-p <remote prefix>                  : Rclone copy to this remote dir prefix (default "")
-r                                  : Remove the file after rclone
-P                                  : Run seeds in parallel
EOF
  exit 1
}

export rclone_remove="false"
export parallel="false"

while getopts c:hde:f:m:o:p:Pr o; do
    case $o in
        (h) usage ;;
        (d) export dryrun="ON";;
        (e) export exec=$OPTARG;;
        (f) export jobfile=$OPTARG;;
        (p) export rclone_prefix=$OPTARG;;
        (r) export rclone_remove="true";;
        (P) export parallel="true";;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done


rclone_file () {
  if [ -z "$rclone_prefix" ]; then
    return
  fi
  if [ -f $1 ] ; then
    if [ -n "$rclone_remove" ] && [ "$rclone_remove" == "true" ]; then
      echo "RCLONE MOVE FILE         : $rclone_prefix $1"
      rclone moveto "$1" neumann:"/mnt/WDB-AN1500/mbl_transition/$rclone_prefix/$1" -L --update
    else
      echo "RCLONE COPY FILE         : $rclone_prefix $1"
      rclone copyto "$1" neumann:"/mnt/WDB-AN1500/mbl_transition/$rclone_prefix/$1" -L --update
    fi
  fi
}

run_sim_id() {
  num_cols=$(awk '{print NF}' $jobfile | head -n 1)
  arg_line=$(tail -n+$1 $jobfile | head -1)
  if [ -z "$arg_line" ]; then
    return 0
  fi
  config_file=$(echo "$arg_line" | cut -d " " -f1)
  config_base=$(echo "$config_file" | xargs -l basename)
  config_dir=$(echo "$config_file" | xargs -l dirname)
  model_seed=$(echo "$arg_line" | cut -d " " -f2)
  outdir=$(awk '$1 ~ /^storage::output_filepath/' $config_file | awk '{sub(/.*=/,""); sub(/ \/!*<.*/,""); print $1;}' | xargs -l dirname)
  outfile=$outdir/mbl_$model_seed.h5
  logdir=logs/$config_dir/$config_base
  logtext=$logdir/$model_seed.txt
  loginfo=$logdir/$model_seed.info
  infoline="SEED:$model_seed|SLURM_ARRAY_JOB_ID:$SLURM_ARRAY_JOB_ID|SLURM_ARRAY_TASK_ID:$SLURM_ARRAY_TASK_ID|SLURM_ARRAY_TASK_STEP:$SLURM_ARRAY_TASK_STEP"
  mkdir -p $logdir

  if [ "$num_cols" -eq 2 ] ; then
    if [ -f $loginfo ] ; then
      echo "Found earlier loginfo: $(tail -n 1 $loginfo)"
      status=$(tail -n 1 $loginfo | awk -F'|' '{print $NF}') # Should be one of RUNNING, FINISHED or FAILED
      if [[ $status =~ FINISHED|RCLONED ]] ; then
        return 0 # Go to next id
      fi
      if [ $status == "RUNNING" ] ; then
        # This could be a simulation that terminated abruptly, or it is actually running right now.
        # We can find out because we can check if the slurm job id is still running using sacct
        old_array_job_id=$(tail -n 1 $loginfo | awk -F'|' '{print $3}' | awk -F':' '{print $2}')
        old_array_task_id=$(tail -n 1 $loginfo | awk -F'|' '{print $4}' | awk -F':' '{print $2}')
        old_job_id=${old_array_job_id}_${old_array_task_id}
        slurm_state=$(sacct -X --jobs $old_job_id --format=state --parsable2 --noheader)
        if [ "$slurm_state" == "RUNNING" ] ; then
          return 0 # Go to next id
        fi
      fi
      # We go a head and run the simulation if it's not running, or if it failed
    fi

    echo "EXEC LINE                : $exec --config=$config_file --outfile=$outfile --seed=$model_seed --threads=$SLURM_CPUS_PER_TASK &>> $logtext"
    if [ -z  "$dryrun" ]; then
      echo "$(date +'%Y-%m-%dT%T')|$infoline|RUNNING" >> $loginfo
      $exec --config=$config_file --outfile=$outfile --seed=$model_seed --threads=$SLURM_CPUS_PER_TASK &>> $logtext
      echo "EXIT CODE                : $?"
      if [ "$?" != "0" ]; then
        echo "$(date +'%Y-%m-%dT%T')|$infoline|FAILED" >> $loginfo
        return $?
      fi
      if [ "$?" == "0" ] ; then
        echo "$(date +'%Y-%m-%dT%T')|$infoline|FINISHED" >> $loginfo
        #logtext='$logdir/$model_seed.txt'
        #outfile=$(awk '/Simulation data written to file/' '$logdir/$model_seed.txt' | awk -F ": " '{print $2}')
        rclone_file $loginfo "false"
        rclone_file $logtext $rclone_remove
        rclone_file $outfile $rclone_remove
        if [ -n "$rclone_prefix" ] && [ "$?" == "0" ]; then
          echo "$(date +'%Y-%m-%dT%T')|$infoline|RCLONED" >> $loginfo
        fi
      fi
    fi
  elif [ "$num_cols" -eq 3 ]; then
    bit_field=$(echo $arg_line | cut -d " " -f3)
    echo "BITFIELD                 : $bit_field"
    echo "EXEC LINE                : $exec -c $config_file -s $model_seed -b $bit_field &>> $logdir/$model_seed_$bit_field.txt"
    if [ -z  "$dryrun" ];then
      $exec -t $SLURM_CPUS_PER_TASK -c $config_file -s $model_seed -b $bit_field &>> $logdir/$model_seed_$bit_field.txt
      exit_code_dmrg=$?
      echo "EXIT CODE         : $exit_code_dmrg"
      return $?
    fi
  else
    echo "Case not implemented"
    exit 1
  fi
}


if [ ! -f $jobfile ]; then
    echo "job file is not a valid file: $jobfile"
    exit 1
fi

echo "HOSTNAME                 : $HOSTNAME"
echo "JOB FILE                 : $jobfile"
echo "SLURM_CLUSTER_NAME       : $SLURM_CLUSTER_NAME"
echo "SLURM_NTASKS             : $SLURM_NTASKS"
echo "SLURM_CPUS_ON_NODE       : $SLURM_CPUS_ON_NODE"
echo "SLURM_JOB_CPUS_PER_NODE  : $SLURM_JOB_CPUS_PER_NODE"
echo "SLURM_CPUS_PER_TASK      : $SLURM_CPUS_PER_TASK" # Task is the same as simulation
echo "SLURM_MEM_PER_CPU        : $SLURM_MEM_PER_CPU"
echo "SLURM_MEM_PER_NODE       : $SLURM_MEM_PER_NODE"
echo "SLURM_JOBID              : $SLURM_JOBID"
echo "SLURM_ARRAY_JOB_ID       : $SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID      : $SLURM_ARRAY_TASK_ID"
echo "SLURM_ARRAY_TASK_STEP    : $SLURM_ARRAY_TASK_STEP"
echo "SLURM_ARRAY_TASK_MIN     : $SLURM_ARRAY_TASK_MIN"
echo "SLURM_ARRAY_TASK_MAX     : $SLURM_ARRAY_TASK_MAX"

ssh-agent -k
eval "$(ssh-agent -s)"

export start_id=$SLURM_ARRAY_TASK_ID
export end_id=$((SLURM_ARRAY_TASK_ID + SLURM_ARRAY_TASK_STEP - 1))
exit_code_save=0
ulimit -c unlimited

echo "TASK ID SEQUENCE         : $(seq -s ' ' $start_id $end_id)"
if [ "$parallel" == "true" ]; then
  # Load GNU Parallel from modules
  module load parallel
  export -f run_sim_id
  export -f rclone_file
  export JOBS_PER_NODE=$SLURM_CPUS_ON_NODE
  if [ -n "$OMP_NUM_THREADS" ]; then
    export JOBS_PER_NODE=$(( $SLURM_CPUS_ON_NODE / $OMP_NUM_THREADS ))
  fi

  parallel --memsuspend=1G --memfree=$SLURM_MEM_PER_CPU \
           --jobs=$JOBS_PER_NODE \
           --ungroup --resume --delay=1s \
           --joblog=logs/parallel-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log \
           --colsep=' ' run_sim_id \
           ::: $(seq $start_id $end_id)
  exit_code_save=$?
else
    for id in $(seq $start_id $end_id); do
      echo "TIME                     : $(date +'%Y-%m-%dT%T')"
      echo "CONFIG LINE              : $arg_line"
      run_sim_id $id
      if [ $? != "0" ]; then
        exit_code_save=$?
      fi
    done
fi

exit $exit_code_save


