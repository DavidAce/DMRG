#!/usr/bin/env bash
rclone delete logs/    --multi-thread-streams=16
rclone delete output/  --multi-thread-streams=16

rclone purge logs/
rclone purge output/

mkdir logs
mkdir output

