/*
 * ABserial.c
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#include "ab_serial.h"
#include "boundBox.h"
#include "computeAuxiliaryGrid.h"
#include "writecell.h"
#include <stdio.h>

void ab_serial(double * xc, double * yc, double * xf, double * yf, int size_c, int size_f){

	double xmin;
	double xmax;
	double ymin;
	double ymax;

	////////////////////////////////////////////////////////////////////
	//		Pre-processing
	////////////////////////////////////////////////////////////////////

	// Create geometry bounding box
	printf("Computing Bounding Box \n");
	boundBox(xf,yf,size_f,&xmin,&xmax,&ymin,&ymax);

	printf("xmin,xmax,ymin,ymax\n");
	printf("%f, %f, %f, %f\n",xmin,xmax,ymin,ymax);


	// Create auxiliary grid
	int resI=10;
	int resJ=10;
	int numCells = (resI-1)*(resJ-1);
	struct cell * auxCells;
	auxCells = (struct cell *)malloc(numCells*sizeof(struct cell));

	computeAuxiliaryGrid(xmin,xmax,ymin,ymax,resI,resJ,auxCells);
	writecell(auxCells,numCells);


}

