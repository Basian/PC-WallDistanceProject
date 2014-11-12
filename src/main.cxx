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
#include "bruteforce/bruteforce.h"

int main(){
	// Grid sizes
	int ni,nj,csize,fsize;

	// cell and face arrays
	double * xc;
	double * yc;
	double * xf;
	double * yf;

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

	// READ GRID AND RETURN CELL AND FACE ARRAYS
	preproc(ni,nj,xc,yc,xf,yf);
	
	// cell and face sizes
	csize = (ni-1)*(nj-1);
	fsize = (ni-1);

	double * wallDist;
	wallDist = new double[csize];

	// BRUTEFORCE SERIAL WALLDISTANCE CALCULATION
	SerialBF(csize,fsize,xc,yc,xf,yf,wallDist);
	
	// WRITE TO FILE
	postproc(ni,nj,wallDist);
	
	// RADIUS CALCULATION
	//double * r;
	//r = new double[csize];
	//radius(csize,xc,yc,r);


	// WRITE TO FILE
	//postproc(ni,nj,r);


	///////////////////////////////////////////////////
	///////////////////////////////////////////////////

	return 0;
}


