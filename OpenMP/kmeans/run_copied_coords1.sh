#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_kmeans

## Output and error files
#PBS -o output-copied-coords1.out
#PBS -e run_kmeans.err

## How many machines should we get? 
#PBS -l nodes=1:ppn=8

##How long should the job run for?
#PBS -l walltime=00:10:00

## Start 
## Run make in the src folder (modify properly)

cd ~/askisi2/kmeans
export GOMP_CPU_AFFINITY="0-63"
for threads in 1 2 4 8 16 32 64
  do
    export OMP_NUM_THREADS=$threads
    echo "Threads: $threads"
    ./kmeans_omp_reduction -s 256 -n 1 -c 4 -l 10 
    echo
  done
