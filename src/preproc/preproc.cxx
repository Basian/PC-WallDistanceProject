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


//void preproc(int &ni, int &nj, double * xc){
void preproc(int &ni, int &nj, double * &xc_buff, double * &yc_buff, double * &xf_buff, double * &yf_buff){
	// grid cell center arrays
	double * xc;
	double * yc;

	// face center arrays
	double * xf;
	double * yf;



	// call fortran subroutine for reading the grid
	// and returning cell/face center x/y arrays
	readmesh_(&ni,&nj,&xc,&yc,&xf,&yf);

	// Allocate buffer arrays for use by MAIN
	xc_buff = new double[(ni-1)*(nj-1)];
	yc_buff = new double[(ni-1)*(nj-1)];

	xf_buff = new double[(ni-1)];
	yf_buff = new double[(ni-1)];


	// Copy data to buffer for use outside of local scope
	for (int i=0; i<(ni-1)*(nj-1); i++){
		xc_buff[i] = xc[i];
		yc_buff[i] = yc[i];
	}

	for (int i=0; i<(ni-1); i++){
		xf_buff[i] = xf[i];
		yf_buff[i] = yf[i];
	}


	return;
}
