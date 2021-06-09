#include<stdio.h>
#include<math.h>
#define F(x) sin(3.1415926*(x)) //初始位移
#define G(x) 0.0 //初始速度
#define a 1.0 //弦长
#define c 2.0 //波动方程系数
#define m 20 //离散时间分割数
#define n 8 //离散空间分割数
#define T 1.0 //模拟时长

int main(int argc,char *argv[]){
    double h; 
    double k; 
    double L; 
    int i,j; //j为时间坐标，i为空间坐标
    double u[m+1][n+1]; //某时刻某点的弦位移
    h=a/n; //空间步
    k=T/m; //时间步
    L=(k*c/h)*(k*c/h); //计算系数
    for (j=0;j<=m;j++) u[j][0]=u[j][n]=0; //两端点位移始终为0
    for (i=1;i<n;i++) u[0][i]=F(i*h); //初始时刻位移
    for (i=1;i<n;i++) //第一次时间迭代后的位移（依赖i-1,i+1,i）
        u[1][i]=(L/2.0)*(u[0][i+1]+u[0][i-1])+(1.0-L)*u[0][i]+k*G(i*h);
    for (j=1;j<m;j++) //时间迭代
        for (i=1;i<n;i++) //空间迭代（依赖j-1,j-2,i-1,i,i+1）
            u[j+1][i]=2.0*(1.0-L)*u[j][i]+L*(u[j][i+1]+u[j][i-1])-u[j-1][i];
    for (j=0;j<=m;j++) //时间扫描
    {
        for (i=0;i<=n;i++)  //空间扫描
            printf("%6.3f\n",u[j][i]); //输出
    }
    return 0;
}
