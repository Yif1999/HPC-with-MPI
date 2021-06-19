
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

//2D矩形平板结构化四边形网格信息
typedef struct  
{
	double dheight;		//矩形高度
	double dwidth;		//矩形宽度	
	double ddelta;		//网格间隔
	int ihei;			//y轴方向网格数目
	int iwid;			//x轴方向网格数目
}SQuadMesh;


//初始化网格数据
void initmeshdata(SQuadMesh mesh, float ***meshval, float **meshval_d, int &irow, int &icol);

void solvethermal(float **meshval, int irow, int icol);

void outputdata(SQuadMesh mesh, float **meshval, int irow, int icol);