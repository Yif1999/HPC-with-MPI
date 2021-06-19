
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

//2D����ƽ��ṹ���ı���������Ϣ
typedef struct  
{
	double dheight;		//���θ߶�
	double dwidth;		//���ο��	
	double ddelta;		//������
	int ihei;			//y�᷽��������Ŀ
	int iwid;			//x�᷽��������Ŀ
}SQuadMesh;


//��ʼ����������
void initmeshdata(SQuadMesh mesh, float ***meshval, float **meshval_d, int &irow, int &icol);

void solvethermal(float **meshval, int irow, int icol);

void outputdata(SQuadMesh mesh, float **meshval, int irow, int icol);