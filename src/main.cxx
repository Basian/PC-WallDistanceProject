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

extern "C" {
#include "advancingBoundary/serial/ab_serial_t2.h"
}

#include "postproc/postproc.h"
#include "bruteforce/bruteforce.h"

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
	
	// cell and face sizes
	size_c = (ni-1)*(nj-1);
	size_f = (ni-1);

	double * wallDistBF;
	wallDistBF = new double[size_c];

	double * wallDistAB;
	wallDistAB = new double[size_c];

	double * wallDistAB_t2;
	wallDistAB_t2 = new double[size_c];

	for(int i=0; i<size_c; i++){
		wallDistAB[i] = 1e9;
		wallDistAB_t2[i] = 1e9;
		wallDistBF[i] = 1e9;
	}

	//////////////////////////////////////////////////
	//////////////////////////////////////////////////

	// BRUTEFORCE SERIAL WALLDISTANCE CALCULATION
	SerialBF(size_c,size_f,xc,yc,xf,yf,wallDistBF);
	
	// WRITE TO FILE
	postproc(ni,nj,wallDistBF,0);
	

	///////////////////////////////////////////////////
	///////////////////////////////////////////////////

	// ABSerial
	ab_serial(xc,yc,xf,yf,size_c,size_f,wallDistAB);
	// WRITE TO FILE
	postproc(ni,nj,wallDistAB,1);


	// ABSerial_t2
	ab_serial_t2(xc,yc,xf,yf,size_c,size_f,wallDistAB_t2);
	// WRITE TO FILE
	postproc(ni,nj,wallDistAB_t2,2);

	///////////////////////////////////////////////////
	///////////////////////////////////////////////////

	return 0;
}


