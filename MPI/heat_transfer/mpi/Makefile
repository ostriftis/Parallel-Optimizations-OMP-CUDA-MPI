CC=mpicc
CFLAGS=-O3 -lm -Wall
#RES=-DPRINT_RESULTS
#CONV=-DTEST_CONV

all: jacobi RedBlack GaussSeidel 

jacobi: Jacobi_mpi.c utils.c
	$(CC) $(CFLAGS) $(RES) $(CONV) Jacobi_mpi.c utils.c -o jacobi_mpi
RedBlack: RedBlackSOR_mpi.c utils.c
	$(CC) $(CFLAGS) $(RES) $(CONV) RedBlackSOR_mpi.c utils.c -o redblackSOR_mpi

GaussSeidel: GaussSeidelSOR_mpi.c utils.c
	$(CC) $(CFLAGS) $(RES) $(CONV) GaussSeidelSOR_mpi.c utils.c -o gaussSeidelSOR_mpi
clean:
	rm jacobi_mpi redblackSOR_mpi gaussSeidelSOR_mpi
