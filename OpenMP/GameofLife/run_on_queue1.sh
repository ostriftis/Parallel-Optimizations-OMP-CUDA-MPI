#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_omp_gameoflife_1

## Output and error files
#PBS -o run_omp_gameoflife_1.out
#PBS -e run_omp_gameoflife_1.err

## How many machines should we get?
#PBS -l nodes=1:ppn=8

##How long should the job run for?
#PBS -l walltime=00:10:00

## Start
## Run make in the src folder (modify properly)

module load openmp
cd /home/parallel/parlab22/askisi1

for array_size in 64 1024 4096
do
  echo "----------------------------------------"
  echo "Array size: $array_size"
  echo "----------------------------------------"
  
  for threads in 1 2 4 6 8
  do
    export OMP_NUM_THREADS=$threads
    echo "Threads: $threads"
    ./omp_gameoflife_1 $array_size 1000
    echo
  done
done
