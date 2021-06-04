#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "MyMPI.h"

typedef double dtype;
#define mpitype MPI_DOUBLE

int main (int argc, char *argv[]) {
  dtype **a; 
  dtype *b;
  dtype *c; 
  dtype *c_part_out; 
  dtype *c_part_in; 
  int *cnt_out;
  int *cnt_in; 
  int *disp_out;
  int *disp_in;
  int i, j;
  int id; 
  int local_els; 
  int m;
  int n;
  int nprime;
  int p;
  dtype *storage; 
  MPI_Init (&argc, &argv) ;
  MPI_Comm_rank (MPI_COMM_WORLD,&id) ;
  MPI_Comm_size (MPI_COMM_WORLD,&p) ;
  read_col_striped_matrix (argv[1], (void***) &a,(void **) &storage, mpitype, &m, &n, MPI_COMM_WORLD) ;
  read_block_vector(argv[2],(void **) &b,mpitype,&nprime,MPI_COMM_WORLD);
  print_block_vector((void *) b,mpitype,nprime,MPI_COMM_WORLD)

  c_part_out=(dtype *) my_malloc (id,n*sizeof(dtype));
  local_els=BLOCK_SIZE(id,p,n);
  for (i=0;i<n;i++){
    c_part_out[i]=0.0;
    for (j=0;j<local_els;j++)
      c_part_out[i]+=a[i][j]*b[j];
  }
  create_mixed_xfer_arrays(id,p,n,&cnt_out,&disp_out);
  create_uniform_xfer_arrays(id,p,n,&cnt_in,&disp_in);
  c_part_in=(dtype *) my_malloc (id,p*local_els*sizeof(dtype));
  MPI_ALLtoallv(c_part_out,cnt_out,disp_out,mpitype,c_part_in,cnt_in,disp_in,mpitype,MPI_COMM_WORLD);
  c=(dtype *) my_malloc (id,local_els*sizeof(dtype));
  for (i=0;i<local_els;i++){
    c[i]=0.0;
    for (j=0;j<p;j++){
      c[i]+=c_part_in[i+j*local_els];
    }
    print_block_vector((void *) c,mpitype,n,MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
  }

