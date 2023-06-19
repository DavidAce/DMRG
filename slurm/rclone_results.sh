#!/usr/bin/env bash
PROGNAME=$0

usage() {
  cat << EOF >&2

Usage              : $PROGNAME [-f <input_file>] [-h] [-m <mode>] [-t <target>]

-a <address>       : Target machine IP address (default = neumann.theophys.kth.se)
-d                 : Performs a dry run.
-h                 : Help. Shows this text.
-p <target prefix> : Prefix at destination (default = /mnt/WDB-AN1500/mbl_transition)
-s <source dir>    : Source relative to current dir (default = .)
-t <target dir>    : Target directory (default = tmp)
-u <user>          : User at target machine (default = david)
-i <include>       : Include this pattern
-L                 : Follow symbolic links
EOF
  exit 1
}
default_adr="neumann.theophys.kth.se"
default_pfx="/mnt/WDB-AN1500/mbl_transition"
default_src="."
default_tgt="tmp"
default_usr="david"
max_threads=8
dry_run=""
follow_sym=""
include=""


while getopts a:d:hi:n:p:s:t:u:L o; do
      case $o in
        (a) default_adr=$OPTARG;;
        (d) dry_run="--dry-run";;
        (h) usage ;;
        (i) include="--include '$OPTARG'";;
        (n) max_threads=$OPTARG;;
        (p) default_pfx=$OPTARG;;
        (s) default_src=$OPTARG;;
        (t) default_tgt=$OPTARG;;
        (u) default_usr=$OPTARG;;
        (L) follow_sym="-L" ;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        (*) usage ;;
      esac
done

#if [ $OPTIND -eq 1 ]; then echo "No flags were passed"; usage ;exit 1; fi


rclone copy $default_src neumann:${default_pfx}/${default_tgt} $follow_sym $dry_run $include --progress --update --multi-thread-streams=$max_threads --transfers=$max_threads
