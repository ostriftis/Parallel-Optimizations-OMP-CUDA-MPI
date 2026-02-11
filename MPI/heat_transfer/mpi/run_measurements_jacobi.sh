#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_mpi_jacobi

## Output and error files
#PBS -o run_mpi_jacobi_measurements_conv.out
#PBS -e run_mpi_jacobi_measurements_conv.err

## How many machines should we get? 
#PBS -l nodes=8:ppn=8


##How long should the job run for?
#PBS -l walltime=00:20:00
 
## Start 
## Run make in the src folder (modify properly)


cd /home/parallel/parlab22/askisi4/heat_transfer/mpi

module load openmpi/1.8.3

proc=(1 2 4 8 16 32 64)
Bx=(1 1 2 2 4 4 8)
By=(1 2 2 4 4 8 8)
sz=(2048 4096 6144)

for s in "${sz[@]}"; do
    
    for ((i=0; i<${#proc[@]}; i++)); do
        
        # Extract values for readability
        p=${proc[i]}
        bx=${Bx[i]}
        by=${By[i]}

        # Calculate Local Sizes
        # Local_X = Global Size / Blocks in X
        # Local_Y = Global Size / Blocks in Y

        # Run MPI
        # Arguments: Global_X Global_Y Local_X Local_Y
        mpirun --mca btl tcp,self -np $p ./jacobi_mpi $s $s $bx $by
    done
done
