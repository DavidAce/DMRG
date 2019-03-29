#!/usr/bin/env bash
mkdir empty_dir
rsync -a -L --delete --progress empty_dir/ logs/
rsync -a -L --delete --progress empty_dir/ output/

rmdir empty_dir