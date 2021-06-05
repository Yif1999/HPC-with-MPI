// 运行时线程数须严格对应维数，单线程处理单行x的迭代
#define TRUE 1
#define FALSE 0
#define bool int

#define MAX_N 100             //允许的最大未知数个数
#define MAX_A (MAX_N * MAX_N) //允许最大的系数的个数

#define MAX_ITERATION 10000 //最大迭代次数
#define TOLERANCE 0.001     //误差

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int pID, pSize;                                //pID：当前进程ID，pSize：总的进程数
int n, iteration = 0;                          //n：未知数的个数，iternation：迭代的次数
float x[MAX_N], new_x[MAX_N], result_x[MAX_N]; //x：表示上一次迭代的结果，new_x：表示这一次迭代的结果，result_x：表示归约之后得到的总的结果
float a[MAX_N][MAX_N];                         //系数
float b[MAX_N];

//输入
void input()
{

    int i, j;

    printf("The n is %d\n", n);
    fflush(stdout);

    //输入系数
    for (i = 0; i < n; i++)
    {
        printf("Input a[%d][0] to a[%d][n-1]:\n", i, i);
        fflush(stdout);
        for (j = 0; j < n; j++)
            scanf("%f", &a[i][j]);
    }

    //输入b
    printf("Input b[0] to b[n-1]:\n");
    fflush(stdout);
    for (j = 0; j < n; j++)
        scanf("%f", &b[j]);
}
//输出结果
void output()
{

    int i;

    printf("Total iteration : %d", iteration);
    if (iteration > MAX_ITERATION)
        printf(", that is out of the limit of iteration!");
    printf("\n");

    for (i = 0; i < n; i++)
        printf("x[%d] is %f\n", i, result_x[i]);
}
//判断精度是否满足要求，满足要求则返回TRUE，否则，返回FALSE
bool tolerance()
{

    int i;

    //有一个不满足误差的，则返回FALSE
    for (i = 0; i < n; i++)
        if (result_x[i] - x[i] > TOLERANCE || x[i] - result_x[i] > TOLERANCE)
            return FALSE;

#ifdef DEBUG
    printf("TRUE From %d\n", pID);
    fflush(stdout);
#endif

    //全部满足误差，返回TRUE
    return TRUE;
}

//入口，主函数
int main(int argc, char *argv[])
{

    MPI_Status status;
    int i, j;
    float sum;

    //初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pID);
    MPI_Comm_size(MPI_COMM_WORLD, &pSize);

    //每个进程对应一个未知数
    n = pSize;

    //根进程负责输入
    if (!pID)
        input();

    //广播系数
    MPI_Bcast(a, MAX_A, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //广播b
    MPI_Bcast(b, MAX_N, MPI_FLOAT, 0, MPI_COMM_WORLD);

#ifdef DEBUG
    //打印a, b
    for (j = 0; j < n; j++)
        printf("%f ", b[j]);
    printf(" From %d\n", pID);
    fflush(stdout);
#endif

    //初始化x的值
    for (i = 0; i < n; i++)
    {
        x[i] = b[i];
        new_x[i] = b[i];
    }

    //分配任务，当前进程处理的元素为该进程的ID
    i = pID;

    //迭代求x[i]的值
    do
    {
        //迭代次数加1
        iteration++;

        //将上一轮迭代的结果作为本轮迭代的起始值，并设置新的迭代结果为0
        for (j = 0; j < n; j++)
        {
            x[j] = result_x[j];
            new_x[j] = 0;
            result_x[j] = 0;
        }

        //同步等待
        MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
        //打印本轮迭代的初始值
        for (j = 0; j < n; j++)
            printf("%f ", x[j]);
        printf(" From %d\n", pID);
        fflush(stdout);
#endif

        // 迭代公式：Dx=-(L+U)x+b，其中D为对角，L为下三角，U为上三角
        //求和，即求(L+U)
        sum = -a[i][i] * x[i];
        for (j = 0; j < n; j++)
            sum += a[i][j] * x[j];

        //求新的x[i]
        new_x[i] = (b[i] - sum) / a[i][i];

#ifdef DEBUG
        //打印本轮迭代的结果
        for (j = 0; j < n; j++)
            printf("%f ", new_x[j]);
        printf(" From %d\n", pID);
        fflush(stdout);
#endif

        //使用归约的方法，同步所有计算结果
        MPI_Allreduce(new_x, result_x, n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

#ifdef DEBUG
        //打印归约后的结果
        for (j = 0; j < n; j++)
            printf("%f ", result_x[j]);
        printf(" From %d\n", pID);
        fflush(stdout);
#endif

        //如果迭代次数超过了最大迭代次数则退出
        if (iteration > MAX_ITERATION)
        {
            break;
        }
    } while (!tolerance()); //精度不满足要求继续迭代

    //根进程负责输出结果
    if (!pID)
        output();

    //结束
    MPI_Finalize();
    return 0;
}