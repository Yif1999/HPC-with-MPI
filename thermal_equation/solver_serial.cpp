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
void initmeshdata(SQuadMesh mesh, float ***meshval, float **meshval_d, int &irow, int &icol)
{
	int i, j;
	double dwidth, dheig;

	irow = mesh.ihei + 1;
	icol = mesh.iwid + 1;

	dwidth = mesh.dwidth;
	dheig = mesh.dheight;

	*meshval_d = new float[irow * icol];

	*meshval = new float *[irow];
	for (i = 0; i < irow; i++)
	{
		(*meshval)[i] = new float[icol];
		for (j = 0; j < icol; j++)
		{
			if (i == irow - 1)
			{
				(*meshval_d)[i * icol + j] = sin(PI * ((float)j / (float)(icol - 1)));
				(*meshval)[i][j] = sin(PI * ((float)j / (float)(icol - 1)));
			}
			else
			{
				(*meshval)[i][j] = 0.0;
				(*meshval_d)[i * icol + j] = 0;
			}
		}
	}
}

#define _WRITE_FILE

/*  �������޲�ַ�����¶ȳ�
 * 	u_new[i,j] = (u_old[i,j-1] + u_old[i,j+1] + u_old[i-1,j] + u_old[i+1,j])/4
 */
void solvethermal(float **meshval, int irow, int icol)
{
	double thre = 1.0; //1.0e-5;
	int i, j, k, iteration = 0;
	float **oldval, eps = 0.0, maxeps = 0.0, global_eps = 0.0;

#ifdef _WRITE_FILE
	FILE *fout = NULL;
	fout = fopen("iter.csv", "w");
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
		printf("Iteration %d  eps:%lf\n", iteration, global_eps);
#ifdef _WRITE_FILE
		fprintf(fout, "%d,%lf\n", iteration, global_eps);
#endif
	} while (global_eps > thre);
	/*while(maxeps > thre);*/

	printf("Convergenced!\n");

#ifdef _WRITE_FILE
	fclose(fout);
#endif

	clock_t t2 = clock();

	double elaps = (double)(t2 - t1) / CLOCKS_PER_SEC;
	printf("CPU elapsed time: %lfs\n", elaps);
}

//����������ݣ�VTK��ʽ�������ڷ���
void outputdata(SQuadMesh mesh, float **meshval, int irow, int icol)
{
	int i, j, iwid, iheig, npts;
	float dwidth, dheig;
	FILE *fout = NULL;

	iheig = mesh.ihei + 1;
	iwid = mesh.iwid + 1;
	npts = iwid * iheig;

	dwidth = mesh.dwidth;
	dheig = mesh.dheight;

	//write VTK file
	fout = fopen("thermal.vtk", "wb");
	if (fout == NULL)
	{
		printf("Error: cann't open file thermal.vtk\n");
		return;
	}
	fprintf(fout, "# vtk DataFile Version 2.0\n");
	fprintf(fout, "2D Thermal file\n");
	fprintf(fout, "ASCII\n");
	fprintf(fout, "DATASET STRUCTURED_GRID\n");
	fprintf(fout, "DIMENSIONS %d %d 1\n", iwid, iheig);
	fprintf(fout, "POINTS %d float\n", npts);
	for (i = 0; i < iheig; i++)
	{
		for (j = 0; j < iwid; j++)
		{
			float x, y;
			x = j * (dwidth / mesh.iwid);
			y = i * (dheig / mesh.ihei);
			fprintf(fout, "%12.6lf %12.6lf 0.0\n", x, y);
		}
	}

	fprintf(fout, "POINT_DATA %d\n", npts);
	fprintf(fout, "SCALARS Thermal_val float\n");
	fprintf(fout, "LOOKUP_TABLE default\n");
	for (i = iheig - 1; i >= 0; i--)
	{
		for (j = 0; j < iwid; j++)
		{
			float val;
			val = meshval[i][j];
			fprintf(fout, "%10.6lf ", val);
			//printf("%10.6lf ", val);
		}
		//printf("\n");
		fprintf(fout, "\n");
	}

	fclose(fout);
	fout = NULL;
}