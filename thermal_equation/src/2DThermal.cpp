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
	MPI_Comm_rank(MPI_COMM_WORLD,&id); //��ȡ���̺�
	MPI_Comm_size(MPI_COMM_WORLD,&p); //��ȡ��������

	if (p==1 && id==0){
		SQuadMesh mesh;
		mesh.dwidth = width;
		mesh.dheight = height;
		mesh.ddelta = delta;	//ͨ����С��ֵ�������������������Ӷ����������ģ(�������ݷֱ�Ϊ0.05��0.015��0.005��0.0015)
		mesh.iwid = mesh.dwidth/mesh.ddelta;
		mesh.ihei = mesh.dheight/mesh.ddelta;

		float **meshval_h, *meshval_d;
		int irow, icol, isize;
		initmeshdata(mesh, &meshval_h, &meshval_d, irow, icol);

		//���м���
		solvethermal_s(meshval_h, irow, icol);
		outputdata(mesh, meshval_h, irow, icol); //��Ҫ���������ʱ�ص�����ע��
	}else if(p>1){
		if (id==0){
			SQuadMesh mesh;
			mesh.dwidth = width;
			mesh.dheight = height;
			mesh.ddelta = delta;	//ͨ����С��ֵ�������������������Ӷ����������ģ(�������ݷֱ�Ϊ0.05��0.015��0.005��0.0015)
			mesh.iwid = mesh.dwidth/mesh.ddelta;
			mesh.ihei = mesh.dheight/mesh.ddelta;

			float **meshval_h, *meshval_d;
			int irow, icol, isize;
			initmeshdata(mesh, &meshval_h, &meshval_d, irow, icol);
		}

		//���м���
		// solvethermal_p(meshval_h, irow, icol);
		// outputdata(mesh, meshval_h, irow, icol);
		printf("emmm...\n");

	}


	
	MPI_Finalize();
	return 0;
}