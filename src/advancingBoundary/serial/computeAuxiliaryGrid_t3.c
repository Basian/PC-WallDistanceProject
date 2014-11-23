/*
 * computeAuxiliaryGrid.c
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#include <stdlib.h>
#include <stdio.h>
#include "computeAuxiliaryGrid_t3.h"

void computeAuxiliaryGrid_t3(double xmin, double xmax, double ymin, double ymax, int resI, int resJ, struct cell_t3 * auxCells){

	int i;
	int j;
	int numCells;

	double * xaux;
	double * yaux;
	xaux = (double *)malloc(resI*sizeof(double));
	yaux = (double *)malloc(resJ*sizeof(double));


	// Compute aux i-array
	for(i=0; i<resI; i++){
		xaux[i] = xmin + ((xmax-xmin)/(double)(resI-1))*(double)i;
	}

	// Compute aux j-array
	for (j=0; j<resJ; j++){
		yaux[j] = ymin + ((ymax-ymin)/(double)(resJ-1))*(double)j;
	}

	// Loop through auxiliary cells and assign local min/max
	for(i=0; i<resI-1; i++){
		for(j=0; j<resJ-1; j++){

			int index = j + i*(resI-1);

			auxCells[index].xmin = xaux[i];
			auxCells[index].xmax = xaux[i+1];
			auxCells[index].ymin = yaux[j];
			auxCells[index].ymax = yaux[j+1];
			auxCells[index].xcenter = (auxCells[index].xmax + auxCells[index].xmin)/2.0;
			auxCells[index].ycenter = (auxCells[index].ymax + auxCells[index].ymin)/2.0;
			auxCells[index].numFaces = 0;
		}
	}


	free(xaux);
	free(yaux);

}
