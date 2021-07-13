#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
  struct eihei
  {
    double x;
    int i;
  };

int main(int argc, char **argv)
{
  int rank,size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

    eihei hi[3][4]={0};

  if (!rank){

    for (int i=0;i<3;i++){
      for (int j=0;j<4;j++){
        hi[i][j].x=j;
      }
    }
    MPI_Send(&hi[2],sizeof(hi)/3/4*3,MPI_BYTE,1,0,MPI_COMM_WORLD);
  }else{
    MPI_Recv(&hi[2],sizeof(hi)/3/4*3,MPI_BYTE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  for (int i=0;i<3;i++){
    for (int j=0;j<4;j++){
          printf("%lf \t",hi[i][j].x);
    }
    printf("\n");
  }

  }

  MPI_Finalize();
 
  return 0;
}
