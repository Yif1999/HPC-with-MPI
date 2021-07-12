#define gama 1.4 //绝热指数
#define Ma 3.0 //来流马赫数
#define C 1.0 //音速
#define P 1.0 //初始压强
#define Rho 1.4 //初始密度
#define TEND 4.0 //解算结束时间
#define dt 0.01 //时间步长
#define dx 0.015//x轴向步长
#define dy 0.015 //y轴向步长
#define  length 3.1 //simBox长度
#define height 1.0 //simBox高度
#define stepL 0.62 //台阶左侧起始位置
#define stepH  0.21 //台阶上部高度位置
#define eps 0.0001 //常数参量

struct velocity
{
    double u;
    double v;
}; //结构体定义速度矢量

struct data 
{
    velocity vel;
    double rho;
    double c;
    double p;
    double E;
    double fp[4];
    double fn[4];
    double gp[4];
    double gn[4];
    double fx[4];
    double gy[4];

}; //结构体定义相关物理量

struct coord
{
   double x;
   double y;
}; //结构体定义坐标

struct unit
{
    int flag=0;
    coord coords;
    data param;
}; //结构体定义有限元

void Block_Divide(int n, int *a); //分块处理函数

unit Steger_Warming_X(unit u);

unit Steger_Warming_Y(unit u);

void WENO_X(double (*fp)[4],double (*fn)[4],double *fx);

void WENO_Y(double (*gp)[4],double (*gn)[4],double *gy);