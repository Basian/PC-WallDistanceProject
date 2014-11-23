/*
 * main.cpp
 *
 *  Created on: Nov 1, 2014
 *      Author: nathan
 */

#include <iostream>
#include <stdio.h>
#include "preproc/preproc.h"
#include "radius/radius.h"

extern "C" {
#include "advancingBoundary/serial/ab_serial.h"
}

extern "C" {
#include "advancingBoundary/serial/ab_serial_t2.h"
}

extern "C" {
#include "advancingBoundary/serial/ab_serial_t3.h"
}

//extern "C" {
#include "advancingBoundary/parallel/ab_pt1.h"
#include "advancingBoundary/parallel/ab_pt2.h"

//}

#include "advancingBoundary/serial/cudaTest.h"


#include "postproc/postproc.h"
#include "bruteforce/bruteforce.h"
#include "bruteforce/bruteforce_parallel.h"


// Utility function for comparing matrices
int compare_matrices(double *out, double *gold, const int N)
{
	int result = 0;
	for (int i = 0; i < N; i++)
	{
		if ( out[i] == gold[i]){}
		else result = 1;
	}
	return result;
}

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

	double * wallDistBFSerial;
	wallDistBFSerial = new double[size_c];

	double * wallDistBFParallel1a;
	wallDistBFParallel1a = new double[size_c];	

	double * wallDistBFParallel1b;
	wallDistBFParallel1b = new double[size_c];	

	double * wallDistBFParallel1c;
	wallDistBFParallel1c = new double[size_c];	

	double * wallDistBFParallel2a;
	wallDistBFParallel2a = new double[size_c];	

	double * wallDistBFParallel2b;
	wallDistBFParallel2b = new double[size_c];	

	double * wallDistBFParallel2c;
	wallDistBFParallel2c = new double[size_c];	

	double * wallDistAB;
	wallDistAB = new double[size_c];

	double * wallDistAB_t2;
	wallDistAB_t2 = new double[size_c];

	double * wallDistAB_t3;
	wallDistAB_t3 = new double[size_c];

	double * wallDistAB_pt1;
	wallDistAB_pt1 = new double[size_c];

	for(int i=0; i<size_c; i++){
		wallDistAB[i] = 1e9;
		wallDistAB_t2[i] = 1e9;
		wallDistAB_t3[i] = 1e9;
		wallDistAB_pt1[i] = 1e9;
		wallDistBFSerial[i] = 1e9;
		wallDistBFParallel1a[i] = 1e9;
		wallDistBFParallel1b[i] = 1e9;
		wallDistBFParallel1c[i] = 1e9;
		wallDistBFParallel2a[i] = 1e9;
		wallDistBFParallel2b[i] = 1e9;
		wallDistBFParallel2c[i] = 1e9;
	}

	//////////////////////////////////////////////////
	//////////////////////////////////////////////////

	// BRUTEFORCE SERIAL WALLDISTANCE CALCULATION
	SerialBF(size_c,size_f,xc,yc,xf,yf,wallDistBFSerial);
	postproc(ni,nj,wallDistBFSerial,0);							// WRITE TO FILE
	
	// BRUTEFORCE PARALLEL 1a WALLDISTANCE CALCULATION
	ParallelBF1a(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel1a);		
	printf("\tVerifying 1a ... %s\n",
	   compare_matrices(wallDistBFParallel1a, wallDistBFSerial, size_c) ? "Failed" : "Success");

   	// BRUTEFORCE PARALLEL 1b WALLDISTANCE CALCULATION
	ParallelBF1b(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel1b);		
	printf("\tVerifying 1b ... %s\n",
	   compare_matrices(wallDistBFParallel1b, wallDistBFSerial, size_c) ? "Failed" : "Success");
		   	
	// BRUTEFORCE PARALLEL 1c WALLDISTANCE CALCULATION
	ParallelBF1c(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel1c);		
	printf("\tVerifying 1c ... %s\n",
	   compare_matrices(wallDistBFParallel1c, wallDistBFSerial, size_c) ? "Failed" : "Success");

   	// BRUTEFORCE PARALLEL 2a WALLDISTANCE CALCULATION
	ParallelBF2a(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel2a);		
	printf("\tVerifying 2a ... %s\n",
	   compare_matrices(wallDistBFParallel2a, wallDistBFSerial, size_c) ? "Failed" : "Success");

   	// BRUTEFORCE PARALLEL 2b WALLDISTANCE CALCULATION
	ParallelBF2b(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel2b);		
	printf("\tVerifying 2b ... %s\n",
	   compare_matrices(wallDistBFParallel2b, wallDistBFSerial, size_c) ? "Failed" : "Success");
	   	
	// BRUTEFORCE PARALLEL 2c WALLDISTANCE CALCULATION
	ParallelBF2c(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel2c);		
	printf("\tVerifying 2c ... %s\n",
	   compare_matrices(wallDistBFParallel2c, wallDistBFSerial, size_c) ? "Failed" : "Success");

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

	// ABSerial_t3
	ab_serial_t3(xc,yc,xf,yf,size_c,size_f,wallDistAB_t3);
	// WRITE TO FILE
	postproc(ni,nj,wallDistAB_t3,2);



	// ABparallel_t1
	ab_parallel_t1(xc,yc,xf,yf,size_c,size_f,wallDistAB_pt1);
	postproc(ni,nj,wallDistAB_pt1,3);

	// ABparallel_t1
	ab_parallel_t2(xc,yc,xf,yf,size_c,size_f,wallDistAB_pt1);
	postproc(ni,nj,wallDistAB_pt1,3);

	///////////////////////////////////////////////////
	///////////////////////////////////////////////////

	return 0;
}


