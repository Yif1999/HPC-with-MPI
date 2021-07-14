#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
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

    /*自适应处理，网格吸附后截断非整除的余数部分*/
    int numL=int(floor(length/dx))/dims[1]*dims[1]; //修正后box长边网格点数
    int numH=int(floor(height/dy))/dims[0]*dims[0]; //修正后box高度网点数
    int numSL=int(floor(stepL/dx)); //修正后台阶前空余网格点数
    int numSH=int(floor(stepH/dy)); //修正后台阶高度网格点数
    float L=(numL-1)*dx; //修正后box长度
    float H=(numH-1)*dy; //修正后box高度
    float sL=(numSL-1)*dx; //修正后台阶前空余长度
    float sH=(numSH-1)*dy; //修正后台阶高度
    
    /*迪卡尔虚拟拓扑*/
    int coords[2];
    MPI_Comm comm2d;
    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&comm2d);
    MPI_Cart_coords(comm2d,id,2,coords); //得到每个进程的迪卡尔坐标存于coords[]

    int blockL=numL/dims[1]+6; //包含虚拟网格的长边网格点数
    int blockH=numH/dims[0]+6; //包含虚拟网格的高度向网格点数
    if (blockL<10||blockH<10){  //若分块中数据过少不足以支撑鬼层则报错
        if (!id)
            printf("[Fatal Error] Too many processes for this problem scale, the execution is abandoned.\n-------------\nTo solve this problem, you can try to reduce the number of  the processes.\n");
        MPI_Finalize();
        return 0;
    }

    /*状态初始化*/
    int i,j;    //行与列的二维索引
    float t;    //记录迭代时间
    float dU[4];    //存储差分导数
    float u_old,v_old,rho_old,E_old,c_old,p_old;    //上一迭代数据缓存
    float corner[2]={sL,sH}; //台阶角点坐标
    unit send_col[(blockH-6)*3],recv_col[(blockH-6)*3]; //通信暂存数据
    //使用malloc在堆上为分块开辟足够大且连续的存储空间
    unit **u=(unit**)malloc(blockH*sizeof(unit*));//分配二维数组
    u[0]=(unit*)malloc(blockH*blockL*sizeof(unit));//分配一维数组
    for(i=1;i<blockH;i++)
        u[i]=u[i-1]+blockL;
    for (i=0;i<blockH;i++){
        #ifdef _OMP_PARALLEL
            #pragma omp parallel for 
        #endif
        for (j=0;j<blockL;j++){
            //基本物理参数赋值
            u[i][j].coords.x=coords[1]*(blockL-6)*dx+(j-3)*dx;
            u[i][j].coords.y=coords[0]*(blockH-6)*dy+(i-3)*dy;
            u[i][j].param.rho=Rho;
            u[i][j].param.vel.u=Ma;
            u[i][j].param.vel.v=0.0;
            u[i][j].param.p=P;
            u[i][j].param.c=C;
            u[i][j].param.E=P/(gama-1)+0.5*Rho*Ma*Ma;
        }
    }

    /*迭代推演*/
    int dest[2];    //目标进程二维坐标
    int rank;   //目标进程号
    int frame=0;    //迭代次数
    float fp[7][4],fn[7][4],gp[7][4],gn[7][4];
    int k,m,n;  //循环索引
    int index;  //一维结果数组索引
    float x_old,y_old;  //边界条件修改前坐标信息暂存
    float res_x[4]={0},res_y[4]={0};    //思维向量暂存差分返回结果
    
    for (t=0;t<=TEND;t+=dt){
        if (!id) printf("%lf percent...\n",t/TEND*100); //打印进度信息

        //边界条件更新
        for (i=0;i<blockH;i++){
            for (j=0;j<blockL;j++){
                //检测模拟框边缘处理边界条件
                if (u[i][j].coords.x>L){
                    x_old=u[i][j].coords.x;
                    u[i][j]=u[i][blockL-4];
                    u[i][j].coords.x=x_old;
                }
                if (u[i][j].coords.y>H){
                    y_old=u[i][j].coords.y;
                    u[i][j]=u[(2*blockH-8-i)][j];
                    u[i][j].param.vel.v=-u[i][j].param.vel.v;
                    u[i][j].coords.y=y_old;
                }
                else if (u[i][j].coords.y<0){
                    y_old=u[i][j].coords.y;
                    u[i][j]=u[6-i][j];
                    u[i][j].param.vel.v=-u[i][j].param.vel.v;
                    u[i][j].coords.y=y_old;
                }
                else if (u[i][j].coords.y==0||u[i][j].coords.y==H){
                    u[i][j].param.vel.v=0.0;
                }
            }
        }

        /*棋盘分块迭代后通信*/
        //Ghost Layer行通信
        if (dims[0]>1){
            if (coords[0]==0){
                dest[0]=coords[0]+1;dest[1]=coords[1];
                MPI_Cart_rank(comm2d,dest,&rank);
                MPI_Sendrecv(&u[blockH-6],sizeof(u)/blockH*3,MPI_BYTE,rank,0,&u[blockH-3],sizeof(u)/blockH*3,MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }else if (coords[0]==dims[0]-1){
                dest[0]=coords[0]-1;dest[1]=coords[1];
                MPI_Cart_rank(comm2d,dest,&rank);
                MPI_Sendrecv(&u[3],sizeof(u)/blockH*3,MPI_BYTE,rank,0,&u[0],sizeof(u)/blockH*3,MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }else{
                dest[0]=coords[0]+1;dest[1]=coords[1];
                MPI_Cart_rank(comm2d,dest,&rank);
                MPI_Sendrecv(&u[blockH-6],sizeof(u)/blockH*3,MPI_BYTE,rank,0,&u[blockH-3],sizeof(u)/blockH*3,MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                dest[0]=coords[0]-1;dest[1]=coords[1];
                MPI_Cart_rank(comm2d,dest,&rank);
                MPI_Sendrecv(&u[3],sizeof(u)/blockH*3,MPI_BYTE,rank,0,&u[0],sizeof(u)/blockH*3,MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        //Ghost Layer列通信
        if (dims[1]>1){
            if (coords[1]==0){
                dest[0]=coords[0];dest[1]=coords[1]+1;
                MPI_Cart_rank(comm2d,dest,&rank);
                for (i=0;i<(blockH-6)*3;i++){
                    send_col[i]=u[3+i/3][blockL-6+i%3];
                }
                MPI_Sendrecv(&send_col,sizeof(send_col),MPI_BYTE,rank,0,&recv_col,sizeof(recv_col),MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for (i=0;i<(blockH-6)*3;i++){
                   u[3+i/3][blockL-3+i%3]= recv_col[i];
                }
            }else if (coords[1]==dims[1]-1){
                dest[0]=coords[0];dest[1]=coords[1]-1;
                MPI_Cart_rank(comm2d,dest,&rank);
                for (i=0;i<(blockH-6)*3;i++){
                    send_col[i]=u[3+i/3][3+i%3];
                }
                MPI_Sendrecv(&send_col,sizeof(send_col),MPI_BYTE,rank,0,&recv_col,sizeof(recv_col),MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);               
                for (i=0;i<(blockH-6)*3;i++){
                   u[3+i/3][i%3]= recv_col[i];
                }
            }else{
                dest[0]=coords[0];dest[1]=coords[1]+1;
                MPI_Cart_rank(comm2d,dest,&rank);
                for (i=0;i<(blockH-6)*3;i++){
                    send_col[i]=u[3+i/3][blockL-6+i%3];
                }
                MPI_Sendrecv(&send_col,sizeof(send_col),MPI_BYTE,rank,0,&recv_col,sizeof(recv_col),MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for (i=0;i<(blockH-6)*3;i++){
                   u[3+i/3][blockL-3+i%3]= recv_col[i];
                }
                dest[0]=coords[0];dest[1]=coords[1]-1;
                MPI_Cart_rank(comm2d,dest,&rank);
                for (i=0;i<(blockH-6)*3;i++){
                    send_col[i]=u[3+i/3][3+i%3];
                }
                MPI_Sendrecv(&send_col,sizeof(send_col),MPI_BYTE,rank,0,&recv_col,sizeof(recv_col),MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);            
                for (i=0;i<(blockH-6)*3;i++){
                   u[3+i/3][i%3]= recv_col[i];
                }
            }
        }

        /*X向角点更新*/
        for (i=0;i<blockH;i++){
            for (j=0;j<blockL;j++){
                if (u[i][j].coords.x>corner[0] && u[i][j].coords.y<corner[1] 
                    && u[i][j].coords.x<=corner[0]+3*dx){   //处于台阶内部位于x轴向虚拟网格内
                    x_old=u[i][j].coords.x;
                    u[i][j]=u[i][j-int(2*(u[i][j].coords.x-corner[0])/dx)];
                    u[i][j].param.vel.u=-u[i][j].param.vel.u;
                    u[i][j].coords.x=x_old;
                }else if (u[i][j].coords.x==corner[0] && u[i][j].coords.y<corner[1]){   //处于台阶边界
                    u[i][j].param.vel.u=0.0;
                }
            }
        }

        /*X向SW通量分解*/
        for (i=0;i<blockH;i++)
            for (j=0;j<blockL;j++)
                u[i][j]=Steger_Warming_X(u[i][j]);

        /*Y向角点更新*/
        for (i=3;i<blockH-3;i++){
            for (j=3;j<blockL-3;j++){
                if (u[i][j].coords.x>corner[0] && u[i][j].coords.y<corner[1] 
                    && u[i][j].coords.y>=corner[1]-3*dy){   //处于台阶内部位于y轴向虚拟网格内
                    y_old=u[i][j].coords.y;
                    u[i][j]=u[i+int(2*(corner[1]-u[i][j].coords.y)/dy)][j];
                    u[i][j].param.vel.v=-u[i][j].param.vel.v;
                    u[i][j].coords.y=y_old;
                }else if(u[i][j].coords.x>corner[0] && u[i][j].coords.y==corner[1]){    //处于台阶边界
                    u[i][j].param.vel.v=0.0;
                }
            }
        }

        /*Y向SW通量分解*/
        for (i=0;i<blockH;i++)
            for (j=0;j<blockL;j++)
                u[i][j]=Steger_Warming_Y(u[i][j]);

        /*棋盘分块角点处理后通信*/
        //Ghost Layer行通信
        if (dims[0]>1){
            if (coords[0]==0){
                dest[0]=coords[0]+1;dest[1]=coords[1];
                MPI_Cart_rank(comm2d,dest,&rank);
                MPI_Sendrecv(&u[blockH-6],sizeof(u)/blockH*3,MPI_BYTE,rank,0,&u[blockH-3],sizeof(u)/blockH*3,MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }else if (coords[0]==dims[0]-1){
                dest[0]=coords[0]-1;dest[1]=coords[1];
                MPI_Cart_rank(comm2d,dest,&rank);
                MPI_Sendrecv(&u[3],sizeof(u)/blockH*3,MPI_BYTE,rank,0,&u[0],sizeof(u)/blockH*3,MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }else{
                dest[0]=coords[0]+1;dest[1]=coords[1];
                MPI_Cart_rank(comm2d,dest,&rank);
                MPI_Sendrecv(&u[blockH-6],sizeof(u)/blockH*3,MPI_BYTE,rank,0,&u[blockH-3],sizeof(u)/blockH*3,MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                dest[0]=coords[0]-1;dest[1]=coords[1];
                MPI_Cart_rank(comm2d,dest,&rank);
                MPI_Sendrecv(&u[3],sizeof(u)/blockH*3,MPI_BYTE,rank,0,&u[0],sizeof(u)/blockH*3,MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        //Ghost Layer列通信
        if (dims[1]>1){
            if (coords[1]==0){
                dest[0]=coords[0];dest[1]=coords[1]+1;
                MPI_Cart_rank(comm2d,dest,&rank);
                for (i=0;i<(blockH-6)*3;i++){
                    send_col[i]=u[3+i/3][blockL-6+i%3];
                }
                MPI_Sendrecv(&send_col,sizeof(send_col),MPI_BYTE,rank,0,&recv_col,sizeof(recv_col),MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for (i=0;i<(blockH-6)*3;i++){
                   u[3+i/3][blockL-3+i%3]= recv_col[i];
                }
            }else if (coords[1]==dims[1]-1){
                dest[0]=coords[0];dest[1]=coords[1]-1;
                MPI_Cart_rank(comm2d,dest,&rank);
                for (i=0;i<(blockH-6)*3;i++){
                    send_col[i]=u[3+i/3][3+i%3];
                }
                MPI_Sendrecv(&send_col,sizeof(send_col),MPI_BYTE,rank,0,&recv_col,sizeof(recv_col),MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);               
                for (i=0;i<(blockH-6)*3;i++){
                   u[3+i/3][i%3]= recv_col[i];
                }
            }else{
                dest[0]=coords[0];dest[1]=coords[1]+1;
                MPI_Cart_rank(comm2d,dest,&rank);
                for (i=0;i<(blockH-6)*3;i++){
                    send_col[i]=u[3+i/3][blockL-6+i%3];
                }
                MPI_Sendrecv(&send_col,sizeof(send_col),MPI_BYTE,rank,0,&recv_col,sizeof(recv_col),MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for (i=0;i<(blockH-6)*3;i++){
                   u[3+i/3][blockL-3+i%3]= recv_col[i];
                }
                dest[0]=coords[0];dest[1]=coords[1]-1;
                MPI_Cart_rank(comm2d,dest,&rank);
                for (i=0;i<(blockH-6)*3;i++){
                    send_col[i]=u[3+i/3][3+i%3];
                }
                MPI_Sendrecv(&send_col,sizeof(send_col),MPI_BYTE,rank,0,&recv_col,sizeof(recv_col),MPI_BYTE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);            
                for (i=0;i<(blockH-6)*3;i++){
                   u[3+i/3][i%3]= recv_col[i];
                }
            }
        }

        /*WENO差分与时间推进*/
        for (i=3;i<blockH-3;i++)
            for (j=3;j<blockL-3;j++){
                //收集单元格差分所需的依赖项信息用于送入WENO函数
                for (m=0;m<7;m++){
                    for (n=0;n<4;n++){
                        fp[m][n]=u[i][j-3+m].param.fp[n];
                        fn[m][n]=u[i][j-3+m].param.fn[n];
                        gp[m][n]=u[i-3+m][j].param.gp[n];
                        gn[m][n]=u[i-3+m][j].param.gn[n];
                    }
                }

                //二维WENO差分计算
                WENO_X(fp,fn,res_x);
                WENO_Y(gp,gn,res_y);

                //计算导数Ut
                for (k=0;k<4;k++){
                    u[i][j].param.fx[k]=res_x[k];
                    u[i][j].param.gy[k]=res_y[k];
                    dU[k]=-(u[i][j].param.fx[k]+u[i][j].param.gy[k]);
                    res_x[k]=0;
                    res_y[k]=0;
                }

                //参数迭代更新
                rho_old=u[i][j].param.rho;
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
                u[i][j].param.c=sqrt(abs(u[i][j].param.p*gama/u[i][j].param.rho));
            }

            //绘制台阶并实时检测是否出现不收敛情形
            for (i=3;i<blockH-3;i++){
                for (j=3;j<blockL-3;j++){
                    if (u[i][j].coords.x>corner[0] && u[i][j].coords.y<corner[1] ){
                        u[i][j].param.rho=-2;
                    }else if (u[i][j].coords.x<corner[0]-dx && u[i][j].coords.y>corner[1]+dy){
                        if (isnan(u[i][j].param.rho)) {
                            printf("[Value Error] Numerical non convergence.\n-------------\nReduce the space and time step size may avoid this problem.\n");
                            MPI_Finalize();
                            return 0;
                        }
                    }
                }
            }

        /*根进程收集数据并导出*/
        if (frame%interval==0){
            //每个进程对所需数据封包
            result perRes[blockH-6][blockL-6];
            for (i=0;i<blockH-6;i++){
                #ifdef _OMP_PARALLEL
                    #pragma omp parallel for 
                #endif
                for (j=0;j<blockL-6;j++){
                    perRes[i][j].x=u[i+3][j+3].coords.x;
                    perRes[i][j].y=u[i+3][j+3].coords.y;
                    perRes[i][j].rho=u[i+3][j+3].param.rho;
                }
            }

            //Gather到根进程
            result totalRes[numH*numL];
            MPI_Gather(perRes,sizeof(perRes),MPI_BYTE,totalRes,sizeof(perRes),MPI_BYTE,0,MPI_COMM_WORLD);

            //写文件
            if (!id){
                char number[5];
                sprintf(number,"%d", frame/interval);
                FILE *fout=NULL;
                char filename[50]={0};
                strcpy(filename,"../out/render");
                strcat(filename,number);
                strcat(filename,".vtk");
                fout=fopen(filename,"wb");
                if (fout == NULL){
                    printf("[Error] cann't open vtk file\n");
                    MPI_Finalize();
                    return 0;
                }
                fprintf(fout, "# vtk DataFile Version 2.0\n");
                fprintf(fout, "2D Sod\n");
                fprintf(fout, "ASCII\n");
                fprintf(fout, "DATASET STRUCTURED_GRID\n");
                fprintf(fout, "DIMENSIONS %d %d 1\n", numL, numH);
                fprintf(fout, "POINTS %d float\n", numL*numH);
                for (i=0;i<numH;i++){
                    for (j=0;j<numL;j++){
                        index=Find_Adrs(i,j,blockL-6,blockH-6,dims);
                        fprintf(fout,"%12.6lf %12.6lf %12.6lf\n",totalRes[index].x,totalRes[index].y,0.0);
                    }
                }
                fprintf(fout, "POINT_DATA %d\n", numL*numH);
                fprintf(fout, "SCALARS rho float\n");
                fprintf(fout, "LOOKUP_TABLE default\n");
                for (i=0;i<numH;i++){
                    for (j=0;j<numL;j++){
                        index=Find_Adrs(i,j,blockL-6,blockH-6,dims);
                        fprintf(fout,"%12.6lf",totalRes[index].rho);
                    }
                    fprintf(fout,"\n");
                }
                fclose(fout);
                fout=NULL;
            }
        }
        frame++;
    }

    //释放内存结束计算
    free(u[0]);
    free(u);
    printf("Finished!\n");
    MPI_Finalize();
    return 0;
}