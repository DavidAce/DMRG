#!/usr/bin/env bash

dir=$1
for i in $(find -L $dir  -maxdepth 1 -type d) ; do
      echo -n $i": " ;
      (find $i  -type f -name *.h5 | wc -l) ;
done

