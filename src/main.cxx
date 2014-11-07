/*
 * main.cpp
 *
 *  Created on: Nov 1, 2014
 *      Author: nathan
 */

#include <iostream>
#include "preproc/preproc.h"
#include "radius/radius.h"
#include "postproc/postproc.h"

int main(){
	// Grid sizes
	int ni,nj,size;

	// cell and face arrays
	double * xc;
	double * yc;
	double * xf;
	double * yf;

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

	// READ GRID AND RETURN CELL AND FACE ARRAYS
	preproc(ni,nj,xc,yc,xf,yf);
	size = (ni-1)*(nj-1);

	double * wallDist;
	wallDist = new double[size];


	// RADIUS CALCULATION
	double * r;
	r = new double[size];
	radius(size,xc,yc,r);


	// WRITE TO FILE
	postproc(ni,nj,r);


	///////////////////////////////////////////////////
	///////////////////////////////////////////////////

	return 0;
}


