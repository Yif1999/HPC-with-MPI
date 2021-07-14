#define gama 1.4 //绝热指数
#define Ma 3.0 //来流马赫数
#define C 1.0 //音速
#define P 1.0 //初始压强
#define Rho 1.4 //初始密度
#define TEND (13*dt)//解算结束时间
#define dt 0.005//时间步长
#define dx 0.02//x轴向步长
#define dy 0.02 //y轴向步长
#define length 3.0 //simBox长度
#define height 1.0 //simBox高度
#define stepL 0.6 //台阶左侧起始位置
#define stepH  0.2 //台阶上部高度位置
#define eps 0.000001 //常数参量

struct result
{
    float rho;
    float x;
    float y;
};

struct velocity
{
    float u;
    float v;
}; //结构体定义速度矢量

struct data 
{
    velocity vel;
    float rho;
    float c;
    float p;
    float E;
    float fp[4];
    float fn[4];
    float gp[4];
    float gn[4];
    float fx[4];
    float gy[4];

}; //结构体定义相关物理量

struct coord
{
   float x;
   float y;
}; //结构体定义坐标

struct unit
{
    int debug=0;
    coord coords;
    data param;
}; //结构体定义有限元

void Block_Divide(int n, int *a); //分块处理函数

unit Steger_Warming_X(unit u);

unit Steger_Warming_Y(unit u);

void WENO_X(float (*fp)[4],float (*fn)[4],float *fx);

void WENO_Y(float (*gp)[4],float (*gn)[4],float *gy);

int Find_Adrs(int i,int j,int blockL,int blockH,int *dims);