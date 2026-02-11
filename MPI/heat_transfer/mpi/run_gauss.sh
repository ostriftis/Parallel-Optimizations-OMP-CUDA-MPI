#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_mpi_gaussSeidel

## Output and error files
#PBS -o run_mpi_gaussSeidel.out
#PBS -e run_mpi_gaussSeidel.err

## How many machines should we get? 
#PBS -l nodes=8:ppn=8


##How long should the job run for?
#PBS -l walltime=00:20:00
 
## Start 
## Run make in the src folder (modify properly)

cd /home/parallel/parlab22/askisi4/heat_transfer/mpi
module load openmpi/1.8.3
export OMPI_MCA_orte_base_help_aggregate=0
mpirun --mca btl tcp,self -np 64 ./gaussSeidelSOR_mpi 512 512 8 8 
