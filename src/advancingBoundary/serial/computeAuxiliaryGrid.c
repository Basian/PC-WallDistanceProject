/*
 * computeAuxiliaryGrid.c
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#include <stdlib.h>
#include "computeAuxiliaryGrid.h"
#include "node.h"




void computeAuxiliaryGrid(int xmin, int xmax, int ymin, int ymax, int resI, int resJ){

	int i,j,numNodes;

	double * xaux;
	double * yaux;
	xaux = (double *)malloc(resI*sizeof(double));
	yaux = (double *)malloc(resJ*sizeof(double));

	numNodes = (resI-1)*(resJ-1);
	struct node * auxNodes;
	auxNodes = (struct node *)malloc(numNodes*sizeof(struct node));


	// Compute aux i-array
	for(i=0; i<resI; i++){
		xaux[i] = xmin + ((xmax-xmin)/(double)resI)*(double)i;
	}

	// Compute aux j-array
	for (j=0; j<resJ; j++){
		yaux[i] = ymin + ((ymax-ymin)/(double)resJ)*(double)j;
	}

	// Loop through auxiliary cells and assign local min/max
	for(i=0; i<resI-1; i++){
		for(j=0; j<resJ-1; j++){

			auxNodes[i*j].xmin = xaux[i];
			auxNodes[i*j].xmax = xaux[i+1];
			auxNodes[i*j].ymin = yaux[j];
			auxNodes[i*j].ymax = yaux[j+1];

		}
	}

}
