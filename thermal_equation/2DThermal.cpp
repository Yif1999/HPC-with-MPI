#include <string.h>
#include "meshdef.h"

int main(int argc, char *argv[])
{
	SQuadMesh mesh;
	mesh.dwidth = 8.0;
	mesh.dheight = 4.0;
	mesh.ddelta = 0.05;	//通过减小此值可以增加网格数量，从而增大问题规模(测试数据分别为0.05、0.015、0.005和0.0015)
	mesh.iwid = mesh.dwidth/mesh.ddelta;
	mesh.ihei = mesh.dheight/mesh.ddelta;

	float **meshval_h, *meshval_d;
	int irow, icol, isize;
	initmeshdata(mesh, &meshval_h, &meshval_d, irow, icol);

	if (strcmp(argv[1],"-s")==0){

	//串行计算
	solvethermal(meshval_h, irow, icol);
	outputdata(mesh, meshval_h, irow, icol); //需要输出计算结果时关掉改行注释
	}else if (strcmp(argv[1],"-p")==0){

	//并行计算



	}
	
	return 0;
}