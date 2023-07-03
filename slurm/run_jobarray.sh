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
-c <config file>                    : Path to a config file (.cfg)
-e <executable>                     : Path to executable (default = "")
-f <jobfile>                        : Path to simulation file, two columns formatted as [configfile seed] (default = "")
-o <seed offset>                    : Start seed count from this offset
-p <remote prefix>                  : Rclone copy to this remote dir prefix (default "")
-r                                  : Remove the file after rclone
-P                                  : Run seeds in parallel
-s                                  : Status file directory
EOF
  exit 1
}
export rclone_remote="neumann:/mnt/WDB-AN1500/mbl_transition"
export rclone_remove="false"
export parallel="false"
export seed_offset=0
export status_dir="status"
while getopts c:hde:f:m:o:p:Prs: o; do
    case $o in
        (h) usage ;;
        (d) export dryrun="ON";;
        (c) printf -v config_path $OPTARG;; # Use printf to expand unicode characters
        (e) export exec=$OPTARG;;
        (o) export seed_offset=$OPTARG;;
        (p) export rclone_prefix=$OPTARG;;
        (r) export rclone_remove="true";;
        (P) export parallel="true";;
        (s) export status_dir=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done
export config_path=$config_path

echodate(){
    printf "%(%Y-%m-%dT%H:%M:%S)T:$*\n" -1
}

log(){
  printf -v statusline "%(%Y-%m-%dT%H:%M:%S)T|$1"
  echo "$statusline" >> $2
}

rclone_copy_to_remote () {
  if [ -z "$rclone_prefix" ]; then
    return
  fi
  if [ -f $1 ] ; then
    if [ "$2" == "true" ]; then
      echodate "RCLONE MOVE LOCAL->REMOTE: $1"
      rclone moveto "$1" "$rclone_remote/$rclone_prefix/$1" -L --update
    else
      echodate "RCLONE COPY LOCAL->REMOTE: $1"
      rclone copyto "$1" "$rclone_remote/$rclone_prefix/$1" -L --update
    fi
  fi
}

rclone_copy_from_remote() {
  if [ -z "$rclone_prefix" ]; then
    return 0
  fi
  file_remote="$rclone_remote/$rclone_prefix/$1"
  file_remote_lsf=$(rclone lsf $file_remote)
  rclone_lsf_exit_code=$?
  if [ "$rclone_lsf_exit_code" == "0" ] && [ "$file_remote_lsf" == "$2" ]; then
    echodate "RCLONE COPY REMOTE->LOCAL: $1"
    rclone copyto $file_remote $1 -L --update
  fi
  return 0 # It's fine if this function fails
}

run_sim_id() {
  array_task_plus_step_id=$1
  model_seed="$(( seed_offset + array_task_plus_step_id ))" # zero-indexed id's

  # Check if there is a status file (which was copied to tmp earlier. Return if finished
  if [ -f "$status_path" ]; then
    status="$(fgrep -e "$model_seed" "$status_path" | cut -d '|' -f2)" # Should get one of TIMEOUT,FAILED,MISSING,FINISHED
    if [ -z "$status" ]; then
          echodate "STATUS                   : NULL $model_seed $id"
          return 0
    fi
    echodate "STATUS                   : $status $model_seed $id"
    if [ "$status" == "FINISHED" ]; then
      return 0
    fi
  fi

  config_file="${config_path##*/}"
  config_base="${config_file%.*}"
  config_dir="${config_path%/*}"
  outdir="${output_path%/*}"
  outfile="$outdir/mbl_$model_seed.h5"
  logdir="logs/$config_dir/$config_base"
  logtext="$logdir/$model_seed.txt"
  loginfo="$logdir/$model_seed.info"
  infoline="SLURM_CLUSTER_NAME:$SLURM_CLUSTER_NAME|HOSTNAME:$HOSTNAME|SEED:$model_seed|SLURM_ARRAY_JOB_ID:$SLURM_ARRAY_JOB_ID|SLURM_ARRAY_TASK_ID:$SLURM_ARRAY_TASK_ID|SLURM_ARRAY_TASK_STEP:$SLURM_ARRAY_TASK_STEP"

  mkdir -p $logdir

  # Next, check if the results already exist in the remote
  # If they do, use rclone copyto to copy the remote file to local
  # This command will only copy if the remote file is newer.
  rclone_copy_from_remote "$loginfo" "$model_seed.info"
  if [ -f $loginfo ] ; then
    echodate "LOGINFO                  : $(tail -n 1 $loginfo)"
    infostatus=$(tail -n 2 $loginfo | awk -F'|' '{print $NF}') # Should be one of RUNNING, FINISHED, RCLONED or FAILED. Add -n 2 to read two lines, in case there is a trailing newline
    if [[ $infostatus =~ FINISHED|RCLONED ]]; then
      # Copy results back to remote
      rclone_copy_to_remote $logtext $rclone_remove
      rclone_copy_to_remote $outfile $rclone_remove
      rclone_exit_code=$?
      if [ -n "$rclone_prefix" ] && [ "$rclone_exit_code" == "0" ]; then
        log "$infoline|RCLONED" "$loginfo"
        rclone_copy_to_remote $loginfo $rclone_remove
      fi
      return 0
    fi

    if [ "$infostatus" == "RUNNING" ] ; then
      # This could be a simulation that terminated abruptly, or it is actually running right now.
      # We can find out because we can check if the slurm job id is still running using sacct
      cluster="$(tail -n 1 $loginfo  | xargs -d '|'  -n1 | grep SLURM_CLUSTER_NAME | awk -F ':' '{print $2}')"
      if [ "$cluster" == "$SLURM_CLUSTER_NAME" ];then
        old_array_job_id="$(tail -n 1 $loginfo  | xargs -d '|'  -n1 | grep SLURM_ARRAY_JOB_ID | awk -F ':' '{print $2}')"
        old_array_task_id="$(tail -n 1 $loginfo  | xargs -d '|'  -n1 | grep SLURM_ARRAY_TASK_ID | awk -F ':' '{print $2}')"
        old_job_id=${old_array_job_id}_${old_array_task_id}
        slurm_state=$(sacct -X --jobs $old_job_id --format=state --parsable2 --noheader)
        if [ "$slurm_state" == "RUNNING" ] ; then
          return 0 # Go to next id
        fi
      elif [ ! -z "$cluster" ]; then
        echodate "WARNING: Job $config_path with seed $model_seed is handled by cluster $cluster"
        return 0 # Go to next id because this job is handled by another cluster
      fi
    fi
    # We go ahead and run the simulation if it's not running, or if it failed
  fi

  # Get the latest data to continue from
  if [ "$status" == "TIMEOUT" ] || [ "$infostatus" == "FAILED" ]; then
    rclone_copy_from_remote "$outfile" "mbl_$model_seed.h5"
    rclone_copy_from_remote "$logtext" "$model_seed.txt"
    extra_args="--revive"
  elif [ "$status" == "FAILED" ];then
    extra_args="--replace"
  fi
  extra_args="--revive"

  echodate "EXEC LINE                : $exec --config=$config_path --outfile=$outfile --seed=$model_seed --threads=$SLURM_CPUS_PER_TASK $extra_args &>> $logtext"
  if [ -z  "$dryrun" ]; then
    trap '$(date +'%Y-%m-%dT%T')|$infoline|FAILED" >> $loginfo' SIGINT SIGTERM
    log "$infoline|RUNNING" "$loginfo"
    $exec --config=$config_path --outfile=$outfile --seed=$model_seed --threads=$SLURM_CPUS_PER_TASK $extra_args &>> $logtext
    exit_code=$?
    echodate "EXIT CODE                : $exit_code"
    if [ "$exit_code" != "0" ]; then
      log "$infoline|FAILED" "$loginfo"
      return $?
    fi
    if [ "$exit_code" == "0" ] ; then
      log "$infoline|FINISHED" "$loginfo"
      rclone_copy_to_remote $logtext $rclone_remove
      rclone_copy_to_remote $outfile $rclone_remove
      rclone_exit_code=$?
      if [ -n "$rclone_prefix" ] && [ "$rclone_exit_code" == "0" ]; then
        log "$infoline|RCLONED" "$loginfo"
        rclone_copy_to_remote $loginfo $rclone_remove
      fi
    fi
  fi
}


if [ ! -f $config_path ]; then
    echodate "config file is not valid: $config_path"
    exit 1
fi

echodate "HOSTNAME                 : $HOSTNAME"
echodate "USER                     : $USER"
echodate "CONFIG PATH              : $config_path"
echodate "SLURM_CLUSTER_NAME       : $SLURM_CLUSTER_NAME"
echodate "SLURM_NTASKS             : $SLURM_NTASKS"
echodate "SLURM_CPUS_ON_NODE       : $SLURM_CPUS_ON_NODE"
echodate "SLURM_JOB_CPUS_PER_NODE  : $SLURM_JOB_CPUS_PER_NODE"
echodate "SLURM_CPUS_PER_TASK      : $SLURM_CPUS_PER_TASK" # Task is the same as simulation
echodate "SLURM_MEM_PER_CPU        : $SLURM_MEM_PER_CPU"
echodate "SLURM_MEM_PER_NODE       : $SLURM_MEM_PER_NODE"
echodate "SLURM_JOBID              : $SLURM_JOBID"
echodate "SLURM_ARRAY_JOB_ID       : $SLURM_ARRAY_JOB_ID"
echodate "SLURM_ARRAY_TASK_ID      : $SLURM_ARRAY_TASK_ID"
echodate "SLURM_ARRAY_TASK_STEP    : $SLURM_ARRAY_TASK_STEP"
echodate "SLURM_ARRAY_TASK_MIN     : $SLURM_ARRAY_TASK_MIN"
echodate "SLURM_ARRAY_TASK_MAX     : $SLURM_ARRAY_TASK_MAX"

if [ -z $SLURM_CLUSTER_NAME ]; then
  echodate "This is not a valid slurm job environment. Use this script with sbatch or srun"
  exit 1
fi

ssh-agent -k
eval "$(ssh-agent -s)"

# Zero-indexed  id's
export start_id=$SLURM_ARRAY_TASK_ID
export end_id=$(( SLURM_ARRAY_TASK_ID + SLURM_ARRAY_TASK_STEP - 1))
exit_code_save=0
ulimit -c unlimited


# Determine the output path from the config file
output_path="$(fgrep -e "storage::output_filepath" "$config_path" | awk '{sub(/.*=/,""); sub(/ \/!*<.*/,""); print $1;}')"
export output_path="$output_path"

# Find and copy the status file to tmp
config_base="$(basename -s ".cfg" "$config_path")"
export status_path="$status_dir/$config_base.status"
if [ -f "$status_path" ]; then
    tempdir="/tmp"
    if [ -d "/scratch/local" ];then
      tempdir="/scratch/local"
    elif [ -n "$PDC_TMP" ]; then
       tempdir="$PDC_TMP"
    fi
    status_temp="$tempdir/DMRG.$USER/status/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    status_name="$config_base.status"
    if [ ! -f "$status_temp/$status_name" ]; then
      mkdir -p $status_temp
      cp $status_path $status_temp/
    fi
    export status_path=$status_temp/$status_name
    trap 'rm -rf "$status_temp"' EXIT
fi

echodate "TASK ID SEQUENCE         : $(seq -s ' ' $start_id $end_id)"
if [ "$parallel" == "true" ]; then
  # Load GNU Parallel from modules
  if [[ "$HOSTNAME" =~ tetralith ]] ; then
    module load parallel/20181122-nsc1
  else
    module load parallel
  fi
  if "$?" != "0"; then
    echo "Failed to module load parallel"
    exit 1
  fi

  export -f echodate
  export -f log
  export -f run_sim_id
  export -f rclone_copy_to_remote
  export -f rclone_copy_from_remote

  export JOBS_PER_NODE=$SLURM_CPUS_ON_NODE
  if [ -n "$OMP_NUM_THREADS" ]; then
    export JOBS_PER_NODE=$(( $SLURM_CPUS_ON_NODE / $OMP_NUM_THREADS ))
  fi

  parallel --memfree=$SLURM_MEM_PER_CPU \
           --jobs=$JOBS_PER_NODE \
           --ungroup --resume --delay=.2s \
           --joblog=logs/parallel-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log \
           --colsep=' ' run_sim_id \
           ::: $(seq $start_id $end_id)
  exit_code_save=$?
else
    for id in $(seq $start_id $end_id); do
      run_sim_id $id
      exit_code=$?
      if [ "$exit_code" != "0" ]; then
        exit_code_save=$exit_code
      fi
    done
fi

exit $exit_code_save


