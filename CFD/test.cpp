#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <util.h>
#include "DynamicArray.h"
 
  struct eihei
  {
    double x;
    int i;
  };

int blockH=3,blockL=5;
void testFunc(double (*a)[blockL][4]){
  printf("%lf\n",a[1][1][1]);
}

int main(int argc, char **argv)
{

  double fp[blockH][blockL][4];
  testFunc(fp);
  

  // int rank,size;

  // MPI_Init(&argc, &argv);
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // MPI_Comm_size(MPI_COMM_WORLD, &size);

  // if (!rank){
  //   printf("%d\n",size);
  //   eihei hi;
  //   hi.x= 0.25;
  //   hi.i=1;
  //   MPI_Send(&hi,sizeof(hi),MPI_BYTE,1,0,MPI_COMM_WORLD);
  // }else{
  //   eihei hello;
  //   hello.x=0;
  //   hello.i=0;
  //   MPI_Recv(&hello,sizeof(hello),MPI_BYTE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  //   printf("%lf %d\n",hello.x,hello.i);
  // }


  // MPI_Finalize();
 
  // return 0;
}
