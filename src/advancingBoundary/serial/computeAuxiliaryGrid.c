/*
 * computeAuxiliaryGrid.c
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#include <stdlib.h>
#include <stdio.h>
#include "computeAuxiliaryGrid.h"

void computeAuxiliaryGrid(int xmin, int xmax, int ymin, int ymax, int resI, int resJ, struct cell * auxCells){

	int i,j,numCells;

	double * xaux;
	double * yaux;
	xaux = (double *)malloc(resI*sizeof(double));
	yaux = (double *)malloc(resJ*sizeof(double));


	printf("xmin,xmax,ymin,ymax\n");
	printf("%f, %f, %f, %f\n",xmin,xmax,ymin,ymax);

	// Compute aux i-array
	for(i=0; i<resI; i++){
//		xaux[i] = xmin + ((xmax-xmin)/(double)resI)*(double)i;
		xaux[i] = xmin + ((xmax-xmin)/resI)*i;
		printf("%f\n",xaux[i]);
	}

	// Compute aux j-array
	for (j=0; j<resJ; j++){
		yaux[i] = ymin + ((ymax-ymin)/(double)resJ)*(double)j;
		printf("%f\n",yaux[i]);
	}

	// Loop through auxiliary cells and assign local min/max
	for(i=0; i<resI-1; i++){
		for(j=0; j<resJ-1; j++){

			auxCells[i*j].xmin = xaux[i];
			auxCells[i*j].xmax = xaux[i+1];
			auxCells[i*j].ymin = yaux[j];
			auxCells[i*j].ymax = yaux[j+1];

		}
	}



}
