#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_kmeans

## Output and error files
#PBS -o check_cpu.out
#PBS -e check_cpu.err

## How many machines should we get?
## #PBS -l nodes=1:ppn=64

##How long should the job run for?
#PBS -l walltime=00:10:00

## Start
## Run make in the src folder (modify properly)

cd ~/askisi2/kmeans
lscpu -e

