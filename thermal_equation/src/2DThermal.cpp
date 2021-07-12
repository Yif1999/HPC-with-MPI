#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "meshdef.h"
#include <mpi.h>
#include <time.h>

int main(int argc, char *argv[])
{
	int id;
	int p;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id); 
	MPI_Comm_size(MPI_COMM_WORLD,&p); 

	SQuadMesh mesh,init;
	mesh.dwidth = 10.0;
	mesh.dheight = 5.0;
	mesh.ddelta = 0.015;	//通过减小此值可以增加网格数量，从而增大问题规模(测试数据分别为0.05、0.015、0.005和0.0015)
	mesh.iwid = mesh.dwidth/mesh.ddelta;
	mesh.ihei = mesh.dheight/mesh.ddelta;
	float **meshval_h;
	int irow, icol, isize;
	init=mesh;

	if (p==1 && id==0){
		initmeshdata_s(mesh, &meshval_h, irow, icol);
		solvethermal_s(meshval_h, irow, icol);		//串行计算
		outputdata(mesh, meshval_h, irow, icol); //输出计算结果
	}else if(p<mesh.ihei+1){
		int rest=(mesh.ihei+1)%p;
		int index=0;
		mesh.ihei=(mesh.ihei+1)/p;
		if (id==p-1){
			mesh.ihei+=(rest+1);
			index=mesh.ihei-rest-1+(id-1)*(mesh.ihei-rest-1);
		}
		else if (id==0)
			mesh.ihei+=1;
		else{
			mesh.ihei+=2;
			index=mesh.ihei-2+(id-1)*(mesh.ihei-2);
		}
		initmeshdata_p(mesh, &meshval_h, irow, icol,id,p);
		clock_t t1 = clock();
		solvethermal_p(mesh, meshval_h, irow, icol, id,p,index);		//并行计算
		clock_t t2 = clock();

		//数据收集，屏蔽鬼层
		int i,j;
    int divide_count[p];
		int block_size=0;
    int total_size=0;
		int disp[p];

		if (id==0 || id==p-1)
			block_size=(irow-1)*icol;
		else
			block_size=(irow-2)*icol;
    MPI_Reduce(&block_size,&total_size,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Gather(&block_size,1,MPI_INT,divide_count,1,MPI_INT,0,MPI_COMM_WORLD);

    disp[0]=0;
    for (i=1;i<p;i++)
      disp[i]=disp[i-1]+divide_count[i-1];

		float unit[block_size];
		if (id==0){
			for (i=0;i<block_size;i++)
				unit[i]=meshval_h[i/icol][i%icol];
		}else{
			for (i=0;i<block_size;i++)
				unit[i]=meshval_h[i/icol+1][i%icol];
		}
		float **res_final;
		res_final = new float *[total_size/icol];
		for (i = 0; i < total_size/icol; i++)
			res_final[i] = new float[icol];
		float res[total_size/icol][icol];
    	MPI_Gatherv(unit,block_size,MPI_FLOAT,res,divide_count,disp,MPI_FLOAT,0,MPI_COMM_WORLD);
		for (i=0;i<total_size/icol;i++){
			for (j=0;j<icol;j++){
				res_final[i][j]=res[i][j];
				// printf("%f ",res_final[i][j]);
			}
			// printf("\n");
		}
		if (id==0) outputdata(init, res_final, total_size/icol, icol); //输出计算结果
		double elaps = (double)(t2 - t1) / CLOCKS_PER_SEC;
		double global_elaps;
		MPI_Reduce(&elaps,&global_elaps,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
		if (id==0){
			printf("CPU elapsed time: %lfs\n", global_elaps);

		}
	}else{
		if (!id)
			printf("[Fatal Error] The number of parallels exceeds the nunber of height partitions, the execution is abandoned.\n-------------------------------------------------------\nTo solve this problem, you can try to reduce the number of parallels.\n");
	}

	MPI_Finalize();
	return 0;
}