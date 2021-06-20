#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "meshdef.h"
#include <mpi.h>

#define width 8.0
#define height 5.0
#define delta 0.05

int main(int argc, char *argv[])
{
	int id;
	int p;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id); //获取进程号
	MPI_Comm_size(MPI_COMM_WORLD,&p); //获取进程总数

	if (p==1 && id==0){
		SQuadMesh mesh;
		mesh.dwidth = width;
		mesh.dheight = height;
		mesh.ddelta = delta;	//通过减小此值可以增加网格数量，从而增大问题规模(测试数据分别为0.05、0.015、0.005和0.0015)
		mesh.iwid = mesh.dwidth/mesh.ddelta;
		mesh.ihei = mesh.dheight/mesh.ddelta;

		float **meshval_h, *meshval_d;
		int irow, icol, isize;
		initmeshdata(mesh, &meshval_h, &meshval_d, irow, icol);

		//串行计算
		solvethermal_s(meshval_h, irow, icol);
		outputdata(mesh, meshval_h, irow, icol); //需要输出计算结果时关掉改行注释
	}else if(p>1){
		if (id==0){
			SQuadMesh mesh;
			mesh.dwidth = width;
			mesh.dheight = height;
			mesh.ddelta = delta;	//通过减小此值可以增加网格数量，从而增大问题规模(测试数据分别为0.05、0.015、0.005和0.0015)
			mesh.iwid = mesh.dwidth/mesh.ddelta;
			mesh.ihei = mesh.dheight/mesh.ddelta;

			float **meshval_h, *meshval_d;
			int irow, icol, isize;
			initmeshdata(mesh, &meshval_h, &meshval_d, irow, icol);
		}

		//并行计算
		// solvethermal_p(meshval_h, irow, icol);
		// outputdata(mesh, meshval_h, irow, icol);
		printf("emmm...\n");

	}


	
	MPI_Finalize();
	return 0;
}