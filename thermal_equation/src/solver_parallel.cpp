#include "meshdef.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define DEBUG

/*
 * x=0, x=dwidth, y=dheight: val = 0;
 * y=0: val[i] = sin((pi*x)/dwidth)
 */
void initmeshdata_p(SQuadMesh mesh, float ***meshval, int &irow, int &icol,int id ,int p){
	int i,j;
	double dwidth,dheig;

	irow=mesh.ihei;
	icol=mesh.iwid+1;

	*meshval=new float *[irow];
	for (i=0;i<irow;i++){
		(*meshval)[i]=new float[icol];
		for (j=0;j<icol;j++){
			(*meshval)[i][j]=0.0;
			if (id==p-1 && i==irow-1){
				(*meshval)[i][j]=sin(PI*((float)j/(float)(icol-1)));
			}
		}
	}

#ifdef DEBUG_
	if (id==2){
		for (i=0;i<irow;i++){
			for (j=0;j<icol;j++)
				printf("%f\t",(*meshval)[i][j]);
			printf("\n");
		}
	}
#endif
}

/*  
 * 	u_new[i,j] = (u_old[i,j-1] + u_old[i,j+1] + u_old[i-1,j] + u_old[i+1,j])/4
 */
void solvethermal_p(SQuadMesh mesh, float **meshval, int irow, int icol,int id ,int p,int index)
{
	double thre = 1.0; //1.0e-5;
	int i, j, k, iteration = 0;
	float **oldval, eps = 0.0, maxeps = 0.0, global_eps = 0.0;
	float global_max=0.0;

	oldval = new float *[irow];
	for (i = 0; i < irow; i++)
		oldval[i] = new float[icol];

	do
	{
		maxeps = 0.0;
		global_eps = 0.0;
		for (i = 0; i < irow; i++)
		{
			for (j = 0; j < icol; j++)
				oldval[i][j] = meshval[i][j];
		}

		for (i = 1; i < irow - 1; i++)
		{
			for (j = 1; j < icol - 1; j++)
			{
				meshval[i][j] = (oldval[i][j - 1] + oldval[i][j + 1] + oldval[i - 1][j] + oldval[i + 1][j]) / 4;

				if (meshval[i][j] == 0.0)
					eps = 10;
				else
					eps = fabs(meshval[i][j] - oldval[i][j]) / meshval[i][j];

				if (fabs(eps) > maxeps)
					maxeps = fabs(eps);
			}

			if (maxeps > global_eps)
				global_eps = maxeps;
		}

		//maxeps = 1.0e-7;
		iteration++;

		// printf("Iteration %d  eps:%lf\n", iteration, global_eps);

#ifdef _WRITE_FILE
		fprintf(fout, "%d,%lf\n", iteration, global_eps);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    if (id==0){
			MPI_Send(meshval[mesh.ihei-2],mesh.iwid+1,MPI_FLOAT,id+1,0,MPI_COMM_WORLD);
			MPI_Recv(meshval[mesh.ihei-1],mesh.iwid+1,MPI_FLOAT,id+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }else if (id==p-1){
			MPI_Send(meshval[1],mesh.iwid+1,MPI_FLOAT,id-1,0,MPI_COMM_WORLD);
			MPI_Recv(meshval[0],mesh.iwid+1,MPI_FLOAT,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }else{
			MPI_Send(meshval[mesh.ihei-2],mesh.iwid+1,MPI_FLOAT,id+1,0,MPI_COMM_WORLD);
			MPI_Recv(meshval[mesh.ihei-1],mesh.iwid+1,MPI_FLOAT,id+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Send(meshval[1],mesh.iwid+1,MPI_FLOAT,id-1,0,MPI_COMM_WORLD);
			MPI_Recv(meshval[0],mesh.iwid+1,MPI_FLOAT,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

#ifdef DEBUG_
	if (id==1){
		for (i=0;i<irow;i++){
			for (j=0;j<icol;j++)
				printf("%f\t",meshval[i][j]);
			printf("\n");
		}
	printf("------------\n");
	}
#endif

	MPI_Allreduce(&global_eps,&global_max,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);

	} while (global_max > thre);
	/*while(maxeps > thre);*/

	// printf("The %d Process Convergenced!\n",id);

#ifdef _WRITE_FILE
	fclose(fout);
#endif

}
