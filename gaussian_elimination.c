// 按行分解demo，未做按块分解
// 运行时线程数须严格对应维数，即一个线程只管一行
#include "stdio.h"
#include "mpi.h"
#include "malloc.h"

#define n 4

int NotIn(int id, int *picked);
struct
{
    double value;
    int index;
} local, global;

int main(int argc, char *argv[])
{
    int id, p;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    double a[n][n + 1] = {{2, 1, -5, 1, 8}, {1, -3, 0, -6, 9}, {0, 2, -1, 2, -5}, {1, 4, -7, 6, 0}};

    int i, j;
    int index;

    int *picked;
    picked = (int *)malloc(sizeof(int) * n); //记录已被选中的行
    for (i = 0; i < n; i++)
        picked[i] = -1;

    for (i = 0; i < n; i++)
    {
        // 利用MPI_MAXLOC选主元
        if (NotIn(id, picked)) //判断该行是否已被选中，没有选择则进行下一步
        {
            local.value = a[id][i];
            local.index = id;
        }
        else
        {
            local.value = -100; //假设各个参数最小值不小于-100
            local.index = id;
        }

        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD); // 归约最大的值，并存入到global中
        // printf(" i = %d,id =%d,value = %f,index = %d\n",i,id,global.value,global.index);
        picked[i] = global.index;

        // 主元对应线程进行广播
        if (id == global.index)
        {
            MPI_Bcast(&global.index, 1, MPI_INT, id, MPI_COMM_WORLD);
        }

        MPI_Allgather(&a[id][0], n + 1, MPI_DOUBLE, a, n + 1, MPI_DOUBLE, MPI_COMM_WORLD); //每个线程解决的是对应行的求解，例如：线程号为0的线程仅得到0行的解，但是第1行的改动，0线程没有办法得到，只有1线程自己才知道，所以需要使用MPI_Allgather（）函数进行去收集，并将结果存入到各个线程中，最后各个线程得到a为最新解

        // 并行消元
        if (NotIn(id, picked))
        {
            double temp = 0 - a[id][i] / a[picked[i]][i];
            for (j = i; j < n + 1; j++)
            {
                a[id][j] += a[picked[i]][j] * temp;
            }
        }
    }

    // 具有先后依赖关系，无法并行，单线程解上三角形
    MPI_Gather(&a[id][0], n + 1, MPI_DOUBLE, a, n + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  //收集各行数据到根线程
    if (id == 0)
    {
        // 打印方程组中的A、b
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n + 1; j++)
            {
                printf("%f\t", a[i][j]);
            }
            printf("\n");
        }
        // 回代求解
        double *x;
        x = (double *)malloc(sizeof(double) * n);
        for (i = (n - 1); i >= 0; i--) 
        {
            x[i] = a[picked[i]][n] / a[picked[i]][i];
            printf("x[%d] = %f\n", i, x[i]);
            for (j = 0; j < n; j++)
            {
                a[picked[j]][n] = a[picked[j]][n] - x[i] * a[picked[j]][i];
                a[picked[j]][i] = 0;
            }
        }
    }

    MPI_Finalize();
    return 0;
}

int NotIn(int id, int *picked)
{
    int i;
    for (i = 0; i < n; i++)
    {
        if (id == picked[i])
        {
            return 0;
        }
    }
    return 1;
}