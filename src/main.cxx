/*
 * main.cpp
 *
 *  Created on: Nov 1, 2014
 *      Author: nathan
 */

#include <iostream>
#include "preproc/preproc.h"
#include "radius/radius.h"

extern "C" {
#include "advancingBoundary/serial/ab_serial.h"
}

#include "postproc/postproc.h"

int main(){
	// Grid sizes
	int ni,nj,size_c,size_f;

	// cell and face arrays
	double * xc;
	double * yc;
	double * xf;
	double * yf;

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

	// READ GRID AND RETURN CELL AND FACE ARRAYS
	preproc(ni,nj,xc,yc,xf,yf);
	size_c = (ni-1)*(nj-1);
	size_f = (ni-1);

	double * wallDist;
	wallDist = new double[size_c];


	// RADIUS CALCULATION
	double * r;
	r = new double[size_c];
	radius(size_c,xc,yc,r);

	// ABSerial
	ab_serial(xc,yc,xf,yf,size_c,size_f);


	// WRITE TO FILE
	postproc(ni,nj,r);


	///////////////////////////////////////////////////
	///////////////////////////////////////////////////

	return 0;
}


