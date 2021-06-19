#include <string.h>
#include "meshdef.h"

int main(int argc, char *argv[])
{
	SQuadMesh mesh;
	mesh.dwidth = 8.0;
	mesh.dheight = 4.0;
	mesh.ddelta = 0.05;	//ͨ����С��ֵ�������������������Ӷ����������ģ(�������ݷֱ�Ϊ0.05��0.015��0.005��0.0015)
	mesh.iwid = mesh.dwidth/mesh.ddelta;
	mesh.ihei = mesh.dheight/mesh.ddelta;

	float **meshval_h, *meshval_d;
	int irow, icol, isize;
	initmeshdata(mesh, &meshval_h, &meshval_d, irow, icol);

	if (strcmp(argv[1],"-s")==0){

	//���м���
	solvethermal(meshval_h, irow, icol);
	outputdata(mesh, meshval_h, irow, icol); //��Ҫ���������ʱ�ص�����ע��
	}else if (strcmp(argv[1],"-p")==0){

	//���м���



	}
	
	return 0;
}