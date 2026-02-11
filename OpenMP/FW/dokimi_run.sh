#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_fw

## Output and error files
#PBS -o dokimi.out
#PBS -e dokimi.err

## How many machines should we get? 
#PBS -l nodes=1:ppn=8

##How long should the job run for?
#PBS -l walltime=00:30:00

## Start 
## Run make in the src folder (modify properly)



cd ~/askisi2/FW

for size in 1024 2048 4096
  do
  echo "Serial"
  ./fw $size
  for threads in 1 2 4 8 16 32 64
    do
      
      export OMP_NUM_THREADS=$threads
      case $threads in
    	1)
      	export GOMP_CPU_AFFINITY="0"
      	;;
    	2)
      	export GOMP_CPU_AFFINITY="0 8"
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
      esac
      echo "GOMP_CPU_AFFINITY=$GOMP_CPU_AFFINITY"
      echo "Threads: $threads"
      echo
      ./fw_sr $size 32
      echo
    done
  done
