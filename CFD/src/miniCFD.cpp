#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define Ma 3.0 //来流马赫数
#define P 1.0 //初始压强
#define Rho 1.4 //初始密度
#define TEND 4.0 //解算结束时间
#define dt 1.0 //时间步长
#define dx 0.1 //x轴向步长
#define dy 0.1 //y轴向步长
#define  length 3.1 //simBox长度
#define height 1.0 //simBox高度
#define stepL 0.62 //台阶左侧起始位置
#define stepH  0.21 //台阶上部高度位置

struct velocity
{
    double u;
    double v;
};

struct data {
    velocity vel;
    double roh;
}; //结构体定义相关物理量

struct coord
{
   double x;
   double y;
}; //结构体定义坐标

struct unit{
    int flag=1;
    coord coords;
    data param[2];
}; //结构体定义有限元

void Block_Divide(int n, int *a); //分块处理函数

int main(int argc, char *argv[]){

    int id,p;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    MPI_Comm_size(MPI_COMM_WORLD,&p);

    int ndims=2;
    int dims[2];
    int periods[2]={0,0};
    int reorder=0;
    Block_Divide(p,dims);   //得到每个棋盘分解块的维度存于dims[]
    // printf("%d %d\n",dims[0],dims[1]);

    /*自适应处理，网格吸附后截断非整除的余数部分*/
    int numL=int(floor(length/dx))/dims[0]*dims[0]; //修正后box长边网格数
    int numH=int(floor(height/dy))/dims[1]*dims[1]; //修正后box高度网格数
    int numSL=int(floor(stepL/dx)); //修正后台阶前空余网格数
    int numSH=int(floor(stepH/dy)); //修正后台阶高度网格数
    double L=numL*dx; //修正后box长度
    double H=numH*dy; //修正后box高度
    double sL=numSL*dx; //修正后台阶前空余长度
    double sH=numSH*dy; //修正后台阶高度
    // printf("%lf %lf %lf %lf\n",L,H,sL,sH);
    
    /*迪卡尔虚拟拓扑*/
    int coords[2];
    MPI_Comm comm2d;
    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&comm2d);
    MPI_Cart_coords(comm2d,id,2,coords); //得到每个进程的迪卡尔坐标存于coords[]
    // printf("%d:%d %d \n",id,coords[0],coords[1]);

    int blockL=numL/dims[0]+6;
    int blockH=numH/dims[1]+6;
    if (blockL<9||blockH<9){
        if (!id)
            printf("[Fatal Error] Too many processes for this problem scale, the execution is abandoned.\n-------------\nTo solve this problem, you can try to reduce the number of  the processes.\n");
        MPI_Finalize();
        return 0;
    }

    /*边界条件初始化*/
    unit u[blockL][blockH]; //为分块开辟存储空间
    int i,j;
    double corner[2]={sL,sH}; //台阶角点坐标
    for (i=0;i<blockH;i++){
        for (j=0;j<blockL;j++){
            u[i][j].coords.x=coords[1]*(blockL-6)*dx+0.5*dx+(j-3)*dx;
            u[i][j].coords.y=coords[0]*(blockH-6)*dy+0.5*dy+(i-3)*dy;
            u[i][j].param[0].roh=Rho;
            u[i][j].param[0].vel.u=Ma*340.0;
            u[i][j].param[0].vel.v=0.0;
            if (u[i][j].coords.x>corner[0] && u[i][j].coords.y<corner[1] 
                && u[i][j].coords.x<corner[0]+3*dx){ //处于台阶内部位于x轴向虚拟网格内
                u[i][j].param[1].vel.u=u[i][j].param[0].vel.u;
                u[i][j].param[0].vel.u=-u[i][j+1-int(2*(u[i][j].coords.x+0.5*dx-corner[0])/dx)].param[0].vel.u;
            }
            // if (id==3) printf("%lf\t%lf\t%lf\t%lf\n",u[i][j].coords->x,u[i][j].coords->y, u[i][j].param->roh, u[i][j].param->u);
        }
    }

    MPI_Finalize();
    return 0;
}





// if (u[i][j].coords.y>corner[1]-3*dy){ //位于y轴向虚拟网格内

//     if (u[i][j].flag){  //标识未被修改，即没有重复

//     }else{  //标识被修改，发生重复，写数据到第二层

//     }
// }