#include "meshdef.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*���ڲ������ֵ��ʼ��Ϊ0���߽������ֵ��Ϊ�߽���������ֵ
 * �߽�����Ϊ��x=0, x=dwidth, y=dheight: val = 0;
 *             y=0: val[i] = sin((pi*x)/dwidth)
 */
void initmeshdata_s(SQuadMesh mesh, float ***meshval, int &irow, int &icol)
{
	int i, j;
	double dwidth, dheig;

	irow = mesh.ihei + 1;
	icol = mesh.iwid + 1;

	dwidth = mesh.dwidth;
	dheig = mesh.dheight;

	*meshval = new float *[irow];
	for (i = 0; i < irow; i++)
	{
		(*meshval)[i] = new float[icol];
		for (j = 0; j < icol; j++)
		{
			if (i == irow - 1)
			{
				(*meshval)[i][j] = sin(PI * ((float)j / (float)(icol - 1)));
			}
			else
			{
				(*meshval)[i][j] = 0.0;
			}
		}
	}
}

// #define _WRITE_FILE

/*  �������޲�ַ�����¶ȳ�
 * 	u_new[i,j] = (u_old[i,j-1] + u_old[i,j+1] + u_old[i-1,j] + u_old[i+1,j])/4
 */
void solvethermal_s(float **meshval, int irow, int icol)
{
	double thre = 0.01; //1.0e-5;
	int i, j, k, iteration = 0;
	float **oldval, eps = 0.0, maxeps = 0.0, global_eps = 0.0;

#ifdef _WRITE_FILE
	FILE *fout = NULL;
	fout = fopen("../out/iter.csv", "w");
	if (!fout)
	{
		printf("Error: cann't open file!\n");
	}
#endif

	clock_t t1 = clock();

	oldval = new float *[irow];
	for (i = 0; i < irow; i++)
		oldval[i] = new float[icol];

	do
	{
		maxeps = 0.0;
		global_eps = 0.0;
		for (i = 0; i < irow; i++)
		{
			for (j = 0; j < icol; j++)
				oldval[i][j] = meshval[i][j];
		}

		for (i = 1; i < irow - 1; i++)
		{
			for (j = 1; j < icol - 1; j++)
			{
				meshval[i][j] = (oldval[i][j - 1] + oldval[i][j + 1] + oldval[i - 1][j] + oldval[i + 1][j]) / 4;

				if (meshval[i][j] == 0.0)
					eps = 10;
				else
					eps = fabs(meshval[i][j] - oldval[i][j]) / meshval[i][j];

				if (fabs(eps) > maxeps)
					maxeps = fabs(eps);
			}

			if (maxeps > global_eps)
				global_eps = maxeps;
		}

		//maxeps = 1.0e-7;
		iteration++;
		//���ÿһ�������в�
		// printf("Iteration %d  eps:%lf\n", iteration, global_eps);
#ifdef _WRITE_FILE
		fprintf(fout, "%d,%lf\n", iteration, global_eps);
#endif
	} while (global_eps > thre);
	/*while(maxeps > thre);*/

	// printf("Convergenced!\n");

#ifdef _WRITE_FILE
	fclose(fout);
#endif

	clock_t t2 = clock();

	double elaps = (double)(t2 - t1) / CLOCKS_PER_SEC;
	printf("CPU elapsed time: %lfs\n", elaps);
}
