#!/bin/bash

## Give the Job a descriptive name
#PBS -N make_omp_gameoflife_1

## Output and error files
#PBS -o make_omp_gameoflife_1.out
#PBS -e make_omp_gameoflife_1.err

## How many machines should we get?
#PBS -l nodes=1:ppn=1

##How long should the job run for?
#PBS -l walltime=00:10:00

## Start
## Run make in the src folder (modify properly)

module load openmp
cd /home/parallel/parlab22/askisi1
make
