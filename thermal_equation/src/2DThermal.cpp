#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "meshdef.h"
#include <mpi.h>

#define width 5.0
#define height 9.0
#define delta 0.05

int main(int argc, char *argv[])
{
	int id;
	int p;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id); 
	MPI_Comm_size(MPI_COMM_WORLD,&p); 

	SQuadMesh mesh;
	mesh.dwidth = width;
	mesh.dheight = height;
	mesh.ddelta = delta;	//通过减小此值可以增加网格数量，从而增大问题规模(测试数据分别为0.05、0.015、0.005和0.0015)
	mesh.iwid = mesh.dwidth/mesh.ddelta;
	mesh.ihei = mesh.dheight/mesh.ddelta;
	float **meshval_h;
	int irow, icol, isize;

	if (p==1 && id==0){
		initmeshdata_s(mesh, &meshval_h, irow, icol);
		solvethermal_s(meshval_h, irow, icol);		//串行计算
		outputdata(mesh, meshval_h, irow, icol); //输出计算结果
	}else if(p<mesh.ihei+1){
		initmeshdata_p(mesh, &meshval_h, irow, icol);

		// solvethermal_p(meshval_h, irow, icol);		//串行计算
		// outputdata(mesh, meshval_h, irow, icol); //输出计算结果
		printf("emmm...\n");

	}else{
		if (!id)
			printf("[Fatal Error] The number of parallels exceeds the nunber of height partitions, the execution is abandoned.\n-------------------------------------------------------\nTo solve this problem, you can try to reduce the number of parallels.\n");
	}

	MPI_Finalize();
	return 0;
}