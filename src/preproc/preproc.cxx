/*
 * preproc.cpp
 *
 *  Created on: Nov 1, 2014
 *      Author: nathan
 */

#include "preproc.h"
#include <iostream>

// fortran declaration for reading grid
extern"C"{
	void readmesh_(int * ni, int * nj, double ** xc, double ** yc, double ** xf, double ** yf);
}


void preproc(){

	// grid point sizes
	int ni;
	int nj;

	// grid cell center arrays
	double * xc;
	double * yc;

	// face center arrays
	double * xf;
	double * yf;


	// call fortran subroutine for reading the grid
	// and returning cell/face center x/y arrays
	readmesh_(&ni,&nj,&xc,&yc,&xf,&yf);



	return;
}
