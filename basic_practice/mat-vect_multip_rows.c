#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "MyMPI.h"

typedef double dtype;
#define mpitype MPI_DOUBLE
int main (int argc,char *argv[]){
    dtype **a;
    dtype *b;
    dtype *c_block;
    dtype *c;
    dtype *storage;
    int i,j;
    int id;
    int m;
    int n;
    int nprime;
    int p;
    int rows;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    read_row_striped_matrix(argv[1],(void ***) &a,(void **) &storage,mpitype,&m,&n,MPI_COMM_WORLD);
    rows=BLOCK_SIZE(id,p,m);
    print_row_striped_matrix((void **) a,mpitype,m,n,MPI_COMM_WORLD);
    read_replicated_vector(argv[2],(void **) &b,mpitype,&nprime,MPI_COMM_WORLD);
    print_replicated_vector(b,mpitype,nprime,MPI_COMM_WORLD);
    c_block=(dtype *)malloc(n*sizeof(dtype));
    c=(dtype*)malloc(n*sizeof(dtype));
    for (i=0;i<rows;i++){
        c_block[i]=0.0;
        for (j=0;j<n;j++)
            c_block[i]+=a[i][j]*b[j];
    }
    replicate_block_vector(c_block,n,(void *) c,mpitype,MPI_COMM_WORLD);
    print_replicated_vector(c,mpitype,n,MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}