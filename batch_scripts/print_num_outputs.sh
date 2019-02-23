#!/usr/bin/env bash
#find -L output/ -type d -print0 | while read -d '' -r dir; do
#    files=("$dir"/*)
#    printf "%5d files in directory %s\n" "${#files[@]}" "$dir"
#done

for i in $(find -L output  -maxdepth 1 -type d) ; do
      echo -n $i": " ;
      (find $i -type f | wc -l) ;
done
