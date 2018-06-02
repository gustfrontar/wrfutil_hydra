#include "mpi.h"
#include <stdio.h>

int main(int argc, char **argv){

  int nprocs, myrank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  printf("dummy-mpi.c nprocs=%d, myrank=%d\n", nprocs, myrank);
  MPI_Finalize();
  return 0;
}
