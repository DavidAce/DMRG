#!/bin/bash



PROGNAME=$0
usage() {
  cat << EOF >&2

Usage                               : $PROGNAME [-options] with the following options:
-h                                  : Help. Shows this text.
-d <path>                           : Directory with failed job analysis (default = "")
-J <str>                            : Job name
-E <end date>                       : End date
-S <start date>                     : Start date
-f <path>                           : Path to failed jobs file
-r <path>                           : Path to resumable jobs file
EOF
  exit 1
}

startdate="-S $(date --iso-8601)"
outfile=failed_jobs.txt
resfile=resume.job
outdir=failed

while getopts hE:d:f:J:r:S: o; do
    case $o in
        (h) usage ;;
        (E) enddate="-E $OPTARG";;
        (d) outdir=$OPTARG;;
        (f) outfile=$OPTARG;;
        (J) jobname="--name $OPTARG";;
        (r) resfile=$OPTARG;;
        (S) startdate="-S $OPTARG";;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done


mkdir -p $outdir
export SACCT_FORMAT="jobid%16,jobname%16,user,state%14,maxrss%10,exitcode%10,account%20,cluster%20"
sacct -u $USER  -X $jobname $startdate $enddate | column -t | tr -s " " | egrep 'OUT|FAI|CAN' > $outdir/$outfile
cat $outdir/$outfile


touch $outdir/$resfile
truncate -s 0 $outdir/$resfile
jobids=$(cat $outdir/$outfile | column -t | cut -d ' ' -f1)
if [ -z "$jobids" ] ; then
    echo "Failed to gather job-ids"
    exit 1
fi

for id in $jobids; do
  logfile=$(find logs/ -name *$id.out)
  if [ -z "$logfile" ]; then
    echo "No match for $id.out"
    continue
  fi
  seed=$(cat $logfile | grep SEED | tail -1 |  cut -d ':' -f2 | xargs)
  cfgfile=$(cat $logfile | grep CONFIG | tail -1 | cut -d ':' -f2 | xargs)
  echo "$cfgfile $seed" >> $outdir/$resfile
  echo "$cfgfile $seed"
done

