#include "meshdef.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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
	fout = fopen("../out/thermal.vtk", "wb");
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
		}
		fprintf(fout, "\n");
	}

	fclose(fout);
	fout = NULL;
}