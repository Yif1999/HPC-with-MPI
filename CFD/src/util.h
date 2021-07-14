#define gama 1.4 //绝热指数
#define Ma 3.0 //来流马赫数
#define C 1.0 //音速
#define P 1.0 //初始压强
#define Rho 1.4 //初始密度
#define TEND 4//解算结束时间
#define dt 0.00005//时间步长
#define dx 0.005//x轴向步长
#define dy 0.005 //y轴向步长
#define length 3.0 //simBox长度
#define height 1.0 //simBox高度
#define stepL 0.6 //台阶左侧起始位置
#define stepH  0.2 //台阶上部高度位置
#define eps 0.000001 //常数参量
#define interval 50 //输出文件的迭代次数间隔
// #define _OMP_PARALLEL //openmp并行开关(注释取消omp/去注释开启omp)

struct result
{
    float rho;
    float x;
    float y;
}; //结构体对最终结果封包

struct velocity
{
    float u;
    float v;
}; //结构体定义速度矢量

struct data 
{
    velocity vel;
    float rho;  //密度
    float c;    //音速
    float p;    //压强
    float E;    //能量
    float fp[4];    //x轴向正分裂
    float fn[4];    //x轴向负分裂
    float fx[4];    //x轴向导数
    float gp[4];    //y轴向正分裂
    float gn[4];    //y轴向负分裂
    float gy[4];    //y轴向导数

}; //结构体定义相关参量

struct coord
{
   float x;
   float y;
}; //结构体定义坐标

struct unit
{
    coord coords;
    data param;
}; //结构体定义有限单元

void Block_Divide(int n, int *a); //棋盘分块维度决策

unit Steger_Warming_X(unit u);  //x轴向SW通量分裂

unit Steger_Warming_Y(unit u);  //y轴向SW通量分裂

void WENO_X(float (*fp)[4],float (*fn)[4],float *fx);   //x轴向WENO差分

void WENO_Y(float (*gp)[4],float (*gn)[4],float *gy);   //y轴向WENO差分

int Find_Adrs(int i,int j,int blockL,int blockH,int *dims); //二维数组映射至Gather后的一维内存地址