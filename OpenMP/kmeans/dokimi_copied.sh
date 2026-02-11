#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_kmeans

## Output and error files
#PBS -o dokimi-copied.out
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
    case $threads in
    1)
        export GOMP_CPU_AFFINITY="0"
        ;;
    2)
        export GOMP_CPU_AFFINITY="0 16"
        ;;
    4)
        export GOMP_CPU_AFFINITY="0 8 16 24"
        ;;
    8)
        export GOMP_CPU_AFFINITY="0 1 8 9 16 17 24 25"
        ;;
    16)
        export GOMP_CPU_AFFINITY="0-3 8-11 16-19 24-27"
        ;;
    32)
        export GOMP_CPU_AFFINITY="0-7 8-15 16-23 24-31"
        ;;
    64)
        export GOMP_CPU_AFFINITY="0-15 16-31 32-47 48-63"
        ;;
    esac
    export OMP_NUM_THREADS=$threads
    echo "Threads: $threads"
    ./kmeans_omp_reduction -s 256 -n 16 -c 32 -l 10 
    echo
  done
