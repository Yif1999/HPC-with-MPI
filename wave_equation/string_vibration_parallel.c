#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

#define F(x) sin(3.1415926*(x)) //初始位移
#define G(x) 0.0 //初始速度
#define a 1.0 //弦长
#define c 2.0 //波动方程系数
#define m 20 //离散时间分割数
#define n 8 //离散空间分割数
#define T 1.0 //模拟时长

/***********Custom Block Communication************/
void Block_Comm(double **u,int j,int block_size,int id,int p){
    if (id==0){
        MPI_Send(&u[j][block_size-2],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD);
        MPI_Recv(&u[j][block_size-1],1 ,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }else if (id==p-1){
        MPI_Send(&u[j][1],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD);
        MPI_Recv(&u[j][0],1 ,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }else{
        MPI_Send(&u[j][block_size-2],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD);
        MPI_Recv(&u[j][block_size-1],1 ,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(&u[j][1],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD);
        MPI_Recv(&u[j][0],1 ,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    return ;
}

/***********Check Ghost Layer************/
int Check_GhLy(int *array, int index,int count){
    for (int i =0; i<count; i++)
        if (array[i]==index)
            return 1;
    return 0;
}

int main(int argc,char *argv[]){
    int id;
    int p;

    double h; 
    double k; 
    double L; 
    int i,j; //j为时间坐标，i为空间坐标
    h=a/n; //空间步
    k=T/m; //时间步
    L=(k*c/h)*(k*c/h); //计算系数

    /**************MPI Initialization****************/
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&id); //获取进程号
    MPI_Comm_size(MPI_COMM_WORLD,&p); //获取进程总数
    if (p>n+1){
        if (!id) printf("[Fatal Error] The number of parallels exceeds the nunber of space partitions, the execution is abandoned.\n-------------------------------------------------------\nTo solve this problem, you can try to reduce the number of parallels.\n");
        MPI_Finalize();
        return 0;
    }else if(p==1){
        if (!id) printf("[Fatal Error] The number of parallels cannot be less than 2, the execution is abandoned.\n-------------------------------------------------------\nTo solve this problem, you can try to increase the number of parallels.\n");
        MPI_Finalize();
        return 0;
    }

    /*************Block Divide**************/
    int block_size=(n+1)/p;
    int index=0;
    //对不同ghost_layer的每一block设定size
    if (id==0)
        block_size+=((n+1)%p+1);
    else if (id==p-1){
        block_size++;
        index=((n+1)/p+(n+1)%p)+(id-1)*((n+1)/p)-1;
    }
    else{
        block_size+=2;
        index=((n+1)/p+(n+1)%p)+(id-1)*((n+1)/p)-1;
    }
    // printf("%d:%d\n",id,index);

    /***********Block Initialization*************/
    double **u=malloc(sizeof(int*)*(m+1));
    for (i=0;i<=m;i++){
        u[i]=malloc(sizeof(double)*block_size);
    }

    //初始时刻位移
    for (i=0;i<block_size;i++){
        u[0][i]=F((index+i)*h);
    }
    //两端位移始终为0
    if (id==0)
        for (j=0;j<=m;j++)
            u[j][0]=0;
    else if (id==p-1)
        for (j=0;j<=m;j++)
            u[j][block_size-1]=0;
    
    /************Iteration************/
    //首次一阶迭代与通信
    for (i=1;i<block_size-1;i++)
        u[1][i]=(L/2.0)*(u[0][i+1]+u[0][i-1])+(1.0-L)*u[0][i]+k*G((index+i)*h);
    MPI_Barrier(MPI_COMM_WORLD);
    Block_Comm(u,1,block_size,id,p);
    
    //循环二阶迭代与通信
    for (j=1;j<m;j++){
        for (i=1;i<block_size-1;i++)
            u[j+1][i]=2.0*(1.0-L)*u[j][i]+L*(u[j][i+1]+u[j][i-1])-u[j-1][i];
        MPI_Barrier(MPI_COMM_WORLD);
        Block_Comm(u,j+1,block_size,id,p);
    }

    /*************Reduction*************/
    int total_size=0;
    MPI_Reduce(&block_size,&total_size,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    int divide_count[p];
    MPI_Gather(&block_size,1,MPI_INT,divide_count,1,MPI_INT,0,MPI_COMM_WORLD);
    int disp[p];
    disp[0]=0;
    for (i=1;i<p;i++)
        disp[i]=disp[i-1]+divide_count[i-1];
    double recv[m+1][total_size];
    for (i=0;i<=m;i++){
        MPI_Gatherv(u[i],block_size,MPI_DOUBLE,recv[i],divide_count,disp,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }


    /*************Root Process Print Resault*************/
    if (!id){
        //计算ghost_layer下表编号存入数组
        int gl_count=total_size-n-1;
        int gl[gl_count];
        for (i=0;i<gl_count;i+=2){
            gl[i]=disp[i/2+1]-1;
            gl[i+1]=disp[i/2+1];
        }
        //打印输出
        for (j=0;j<=m;j++){ 
            for (i=0;i<total_size;i++)  
                    if (!Check_GhLy(gl,i,gl_count)) 
                        printf("%6.3f\t",recv[j][i]); 
            printf("\n");
        }
    }

    /***********Free Malloc************/
    for (i=0;i<m+1;i++)
        free(u[i]);
    free(u);

    MPI_Finalize();
    return 0;
}