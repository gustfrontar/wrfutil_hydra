#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef MPI_ARGV_NULL
#define MPI_ARGV_NULL (NULL)
#endif
#define SPAWNPATHLEN 1024
#define SPAWNPATHLENS "1024"

int main(int argc, char **argv){

  int maxspawn = 100000;
  int nprocs, myrank, new_nprocs[maxspawn], nspawn, i;
  char path[maxspawn][SPAWNPATHLEN], wdir[maxspawn][SPAWNPATHLEN], cwd[SPAWNPATHLEN];
  MPI_Comm intercomm;
  MPI_Info info;
  int *ary_ierr;
  FILE *fp;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  printf("nprocs=%d, myrank=%d\n", nprocs, myrank);

  nspawn = 0;
  if(myrank == 0){
    fp = fopen(argv[1], "r");
    while(feof(fp) == 0){
      fscanf(fp, "%d %" SPAWNPATHLENS "s %" SPAWNPATHLENS "s", &(new_nprocs[nspawn]), path[nspawn], wdir[nspawn]);
      printf("%d %s %s\n", new_nprocs[nspawn], path[nspawn], wdir[nspawn]);
      if(wdir[nspawn][0] == 0) break;
      nspawn += 1;
    }
    printf("nspawn = %d\n", nspawn);
    fclose(fp);
    //if(nprocs != new_nprocs[nspawn - 1]){
    //  printf("Error: rank of spawn.c (%d) != rank of the last command in the input file (%d)\n", nprocs, new_nprocs[nspawn - 1]);
    //  MPI_Abort(MPI_COMM_WORLD, 999);
    //}
  }
  MPI_Bcast(&nspawn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(new_nprocs, maxspawn, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(path, SPAWNPATHLEN * maxspawn, MPI_CHARACTER, 0, MPI_COMM_WORLD);
  MPI_Bcast(wdir, SPAWNPATHLEN * maxspawn, MPI_CHARACTER, 0, MPI_COMM_WORLD);
  
  MPI_Info_create(&info);
  //for(i = 0; i < nspawn - 1; i++){
  for(i = 0; i < nspawn; i++){
    printf("Spawn process %d\n", i);
    ary_ierr = malloc(sizeof(int) * new_nprocs[i]);
    getcwd(cwd, SPAWNPATHLEN);
    if(chdir(wdir[i])){
      printf("Error in chdir(%s)\n", wdir[i]);
      MPI_Abort(MPI_COMM_WORLD, 999);
    }
    MPI_Info_set(info, "wdir", wdir[i]);
    MPI_Comm_spawn(path[i], MPI_ARGV_NULL, new_nprocs[i], info, 0, MPI_COMM_WORLD, &intercomm, ary_ierr);
    chdir(cwd);
    free(ary_ierr);
  }
  MPI_Info_free(&info);
  MPI_Finalize();

  /* exec the last process */
  //if(chdir(wdir[nspawn - 1])){
  //  printf("Error in chdir(%s)\n", wdir[nspawn - 1]);
  //  exit(1);
  //}
  //execl(path[nspawn - 1], path[nspawn - 1], NULL);

  /* these may be executed only if execl fails */
  //printf("spawn.c: Error in running the last command in the list.\n");
  return 0;
}
