#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

//2D平板网格结构定义
typedef struct  
{
	double dheight;		//平板高度
	double dwidth;		//平板宽度
	double ddelta;		//细分粒度
	int ihei;			//高度格数
	int iwid;			//宽度格数
}SQuadMesh;


//声明函数
void initmeshdata_s(SQuadMesh mesh, float ***meshval, int &irow, int &icol);

void solvethermal_s(float **meshval, int irow, int icol);

void initmeshdata_p(SQuadMesh mesh, float ***meshval, int &irow, int &icol,int id,int p);

void solvethermal_p(SQuadMesh mesh, float **meshval, int irow, int icol, int id,int p,int index);

void outputdata(SQuadMesh mesh, float **meshval, int irow, int icol);