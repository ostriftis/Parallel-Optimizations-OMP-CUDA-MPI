#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_mpi_kmeans

## Output and error files
#PBS -o run_mpi_kmeans.out
#PBS -e run_mpi_kmeans.err

## How many machines should we get? 
#PBS -l nodes=8:ppn=8


##How long should the job run for?
#PBS -l walltime=00:20:00
 
## Start 
## Run make in the src folder (modify properly)

cd /home/parallel/parlab22/askisi4/kmeans
module load openmpi/1.8.3

for i in 1 2 4 8 16 32 64; do
  mpirun --mca btl tcp,self -np $i ./kmeans_mpi -s 256 -n 16 -c 32 -l 10
done
