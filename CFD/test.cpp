#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
int main(int argc, char **argv)
{
  int rank,size;
 
  int rankX, rankY;
  int ndims = 2;
  int dims[2] = {2, 2}; 
  int periods[2] = {0, 0}; 
  int reorder = 0;
 
  int remainX[2] = {1, 0}; 
  int remainY[2] = {0, 1}; 
 
  MPI_Comm comm2d;
  MPI_Comm commX, commY;
 
  MPI_Init(&argc, &argv);
 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
 
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm2d);
  MPI_Cart_sub(comm2d, remainX, &commX);
  MPI_Cart_sub(comm2d, remainY, &commY);
 
  MPI_Comm_rank(commX, &rankX);
  MPI_Comm_rank(commY, &rankY);
 
  printf("rank = %d;   X = %d;  Y = %d\n", rank, rankX, rankY);
 
  MPI_Finalize();
 
  return 0;
}
