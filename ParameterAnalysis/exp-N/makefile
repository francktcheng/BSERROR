#compile cpu mpi code
par_mpi.cpu : par_mpi.c inc/onetimesimu.h inc/misc.h
	mpiicc -Wall -std=c99 -mkl -openmp par_mpi.c -O2 -o par_mpi.cpu 

