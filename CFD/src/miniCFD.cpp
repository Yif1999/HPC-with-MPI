#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include "util.h"

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

    int blockL=numL/dims[0]+6; //包含虚拟网格的长边网格数
    int blockH=numH/dims[1]+6; //包含虚拟网格的高度向网格数
    if (blockL<9||blockH<9){
        if (!id)
            printf("[Fatal Error] Too many processes for this problem scale, the execution is abandoned.\n-------------\nTo solve this problem, you can try to reduce the number of  the processes.\n");
        MPI_Finalize();
        return 0;
    }

    /*状态初始化*/
    printf("%d %d\n",blockH,blockL);
    unit u[blockH][blockL]; //为分块开辟存储空间
    printf("eihei\n");
    fflush(stdout);
    int i,j;
    double t;
    double dU[4];
    double u_old,v_old,rho_old,E_old,c_old,p_old;
    double corner[2]={sL,sH}; //台阶角点坐标
    for (i=0;i<blockH;i++){
        for (j=0;j<blockL;j++){
            //基本物理参数赋值
            u[i][j].coords.x=coords[1]*(blockL-6)*dx+0.5*dx+(j-3)*dx;
            u[i][j].coords.y=coords[0]*(blockH-6)*dy+0.5*dy+(i-3)*dy;
            u[i][j].param.rho=Rho;
            u[i][j].param.vel.u=Ma;
            u[i][j].param.vel.v=0.0;
            u[i][j].param.p=P;
            u[i][j].param.c=C;
            u[i][j].param.E=P/(gama-1)+0.5*Rho*Ma*Ma;
        }
    }
    /*迭代推演*/
    double fp[7][4],fn[7][4],gp[7][4],gn[7][4];
    int k,m,n;
    double res_x[4],res_y[4];
    for (t=0;t<=4;t+=dt){
        // if (!id) printf("%lf percent...\n",t/4*100);

        //边界条件更新
        for (i=0;i<blockH;i++){
            for (j=0;j<blockL;j++){
                //检测模拟框边缘
                if (u[i][j].coords.x>L){
                    u[i][j]=u[i][blockL-4];
                }
                if (u[i][j].coords.y>H){
                    u[i][j]=u[(2*blockH-7-i)][j];
                    u[i][j].param.vel.v=-u[i][j].param.vel.v;
                }
                else if (u[i][j].coords.y<0){
                    u[i][j]=u[5-i][j];
                    u[i][j].param.vel.v=-u[i][j].param.vel.v;
                }
            }
        }

        /*X向角点更新*/
        for (i=0;i<blockH;i++){
            for (j=0;j<blockL;j++){
                if (u[i][j].coords.x>corner[0] && u[i][j].coords.y<corner[1] 
                    && u[i][j].coords.x<corner[0]+3*dx){   //处于台阶内部位于x轴向虚拟网格内
                    u[i][j]=u[i][j+1-int(2*(u[i][j].coords.x+0.5*dx-corner[0])/dx)];
                    u[i][j].param.vel.u=u[i][j].param.vel.u;
                }
            }
        }

        /*X向SW通量分解*/
        for (i=0;i<blockH;i++)
            for (j=0;j<blockL;j++)
                u[i][j]=Steger_Warming_X(u[i][j]);

        /*Y向角点更新*/
        for (i=0;i<blockH;i++){
            for (j=0;j<blockL;j++){
                if (u[i][j].coords.x>corner[0] && u[i][j].coords.y<corner[1] 
                    && u[i][j].coords.y>corner[1]-3*dy){   //处于台阶内部位于y轴向虚拟网格内
                    u[i][j]=u[i+1-int(2*(u[i][j].coords.y+0.5*dy-corner[1])/dy)][j];
                    u[i][j].param.vel.v=u[i][j].param.vel.v;
                }
            }
        }

        /*Y向SW通量分解*/
        for (i=0;i<blockH;i++)
            for (j=0;j<blockL;j++)
                u[i][j]=Steger_Warming_Y(u[i][j]);

        /*WENO差分与时间推进*/
        for (i=3;i<blockH-3;i++)
            for (j=3;j<blockL-3;j++){
                for (m=0;m<7;m++){
                    for (n=0;n<4;n++){
                        fp[m][n]=u[i][j-3+m].param.fp[n];
                        fn[m][n]=u[i][j-3+m].param.fn[n];
                        gp[m][n]=u[i-3+m][j].param.gp[n];
                        gn[m][n]=u[i-3+m][j].param.gn[n];
                    }
                }

                WENO_X(fp,fn,res_x);
                WENO_Y(gp,gn,res_y);
                for (k=0;k<4;k++){
                    u[i][j].param.fx[k]=res_x[k];
                    u[i][j].param.gy[k]=res_y[k];
                    dU[k]=-(u[i][j].param.fx[k]+u[i][j].param.gy[k]);
                }

                //参数迭代更新
                rho_old=u[i][j].param.rho;
                printf("%lf\n",rho_old);
                E_old=u[i][j].param.E;
                c_old=u[i][j].param.c;
                p_old=u[i][j].param.p;
                u_old=u[i][j].param.vel.u;
                v_old=u[i][j].param.vel.v;
                u[i][j].param.rho=rho_old+dt*dU[0];
                u[i][j].param.vel.u=(rho_old*u_old+dt*dU[1])/u[i][j].param.rho;
                u[i][j].param.vel.v=(rho_old*v_old+dt*dU[2])/u[i][j].param.rho;
                u[i][j].param.E=E_old+dt*dU[3];
                u[i][j].param.p=(gama-1)*(u[i][j].param.E-0.5*u[i][j].param.rho*(u[i][j].param.vel.u*u[i][j].param.vel.u+u[i][j].param.vel.v*u[i][j].param.vel.v));
                u[i][j].param.c=sqrt(u[i][j].param.p*gama/u[i][j].param.rho);
            }
    }

    if (0){
        FILE *fout=NULL;
        fout=fopen("../out/test.vtk","wb");
        if (fout == NULL){
            printf("Error: cann't open file thermal.vtk\n");
            MPI_Finalize();
            return 0;
        }
        fprintf(fout, "# vtk DataFile Version 2.0\n");
        fprintf(fout, "2D Thermal file\n");
        fprintf(fout, "ASCII\n");
        fprintf(fout, "DATASET STRUCTURED_GRID\n");
        fprintf(fout, "DIMENSIONS %d %d 1\n", numL, numH);
        fprintf(fout, "POINTS %d float\n", numL*numH);
        for (i=3;i<blockH-3;i++){
            for (j=3;j<blockL-3;j++){
                fprintf(fout,"%12.6lf %12.6lf %12.6lf\n",u[i][j].coords.x,u[i][j].coords.y,0.0);
            }
        }
        fprintf(fout, "POINT_DATA %d\n", numL*numH);
        fprintf(fout, "SCALARS Thermal_val float\n");
        fprintf(fout, "LOOKUP_TABLE default\n");
        for (i=3;i<blockH-3;i++){
            for (j=3;j<blockL-3;j++){
                fprintf(fout,"%12.6lf",u[i][j].param.rho);
            }
            fprintf(fout,"\n");
        }
        fclose(fout);
        fout=NULL;
    }




    MPI_Finalize();
    return 0;
}

