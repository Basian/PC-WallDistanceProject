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

#include "advancingBoundary/serial/cudaTest.h"


#include "postproc/postproc.h"
#include "bruteforce/bruteforce.h"
#include "bruteforce/bruteforce_parallel.h"
#include <sys/resource.h>
#include <sys/times.h>

// Utility function for computing elapsed time
double gettime_msec(struct rusage *tm_start, struct rusage *tm_end)
{
	double end = (double)tm_end->ru_utime.tv_sec + (double)tm_end->ru_utime.tv_usec / 1000000.0;
	double start = (double)tm_start->ru_utime.tv_sec + (double)tm_start->ru_utime.tv_usec / 1000000.0;

	double diff = end-start;
	
	return diff*1000; // return time in msec

}

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

	// define timer variables for calculating time
	struct rusage tm_start;
	struct rusage tm_end;
	double time_in_msec;

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

	// READ GRID AND RETURN CELL AND FACE ARRAYS
	preproc(ni,nj,xc,yc,xf,yf);
	
	// cell and face sizes
	size_c = (ni-1)*(nj-1);
	size_f = (ni-1);

	double * wallDistBFSerial;
	wallDistBFSerial = new double[size_c];

	double * wallDistBFParallel1;
	wallDistBFParallel1 = new double[size_c];	

	double * wallDistBFParallel2;
	wallDistBFParallel2 = new double[size_c];	

	double * wallDistAB;
	wallDistAB = new double[size_c];

	double * wallDistAB_t2;
	wallDistAB_t2 = new double[size_c];

	for(int i=0; i<size_c; i++){
		wallDistAB[i] = 1e9;
		wallDistAB_t2[i] = 1e9;
		wallDistBFSerial[i] = 1e9;
		wallDistBFParallel1[i] = 1e9;
		wallDistBFParallel2[i] = 1e9;
	}

	//////////////////////////////////////////////////
	//////////////////////////////////////////////////

	// BRUTEFORCE SERIAL WALLDISTANCE CALCULATION
	getrusage( RUSAGE_SELF, &tm_start ); 						// Start timer
	SerialBF(size_c,size_f,xc,yc,xf,yf,wallDistBFSerial);
	getrusage( RUSAGE_SELF, &tm_end );   						// End timer
	time_in_msec = gettime_msec( &tm_start, &tm_end ); 			// Get elapsed time
	
	printf("Brute force - serial execution time: %g milliseconds.\n", time_in_msec);
	postproc(ni,nj,wallDistBFSerial,0);							// WRITE TO FILE
	
	// BRUTEFORCE PARALLEL 1 WALLDISTANCE CALCULATION
	getrusage( RUSAGE_SELF, &tm_start ); 						// Start timer
	ParallelBF1(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel1);		
	getrusage( RUSAGE_SELF, &tm_end );   						// End timer
	time_in_msec = gettime_msec( &tm_start, &tm_end ); 			// Get elapsed time

	printf("Brute force - parallel 1 execution time (block per cell): %g milliseconds.\n\tVerifying output result ...%s\n",
		   time_in_msec, compare_matrices(wallDistBFParallel1, wallDistBFSerial, size_c) ? "Failed" : "Success");

   	// BRUTEFORCE PARALLEL 2 WALLDISTANCE CALCULATION
	getrusage( RUSAGE_SELF, &tm_start ); 						// Start timer
	ParallelBF2(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel2);		
	getrusage( RUSAGE_SELF, &tm_end );   						// End timer
	time_in_msec = gettime_msec( &tm_start, &tm_end ); 			// Get elapsed time

	printf("Brute force - parallel 2 execution time (block per cell shared mem): %g milliseconds.\n\tVerifying output result ...%s\n",
		   time_in_msec, compare_matrices(wallDistBFParallel2, wallDistBFSerial, size_c) ? "Failed" : "Success");
		   	
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
	cudaTest();


	return 0;
}


