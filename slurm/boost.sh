#!/bin/bash

for job in $(squeue -u x_davac | cut -d " " -f1); do
        nsc-boost-priority $job
done
