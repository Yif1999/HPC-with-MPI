#include <mpi.h>
#include <stdio.h>
#include <omp.h>

int main(int argc,char *argv[]){
  int i;
  int id;
  int p;
  double count;
  double global_count;
  int n;
  double one_step(int,int);

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  MPI_Comm_size(MPI_COMM_WORLD,&p);

  n=10000000;count=0;global_count=0;

  #pragma omp parallel for
  for (i=id;i<n;i+=p){
    double d=one_step(n,i);
    #pragma omp critical
    count+=d;
  }
  
  // printf("Process %d is done\n",id);
  // fflush(stdout);
  MPI_Reduce (&count,
              &global_count,
              1,
              MPI_DOUBLE,
              MPI_SUM,
              0,
              MPI_COMM_WORLD);
  if (id==0) printf("pi=%lf\n",global_count);
  MPI_Finalize();

  return 0;  
}

double one_step(int n,int i){
  double x;
  double area;
  x=i*1.0/n;
  area=1.0/n*(4/(1+x*x));
  return area;
}