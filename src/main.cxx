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

//extern "C" {
#include "advancingBoundary/parallel/ab_pt1.h"
//}

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

	double * wallDistBFParallel3;
	wallDistBFParallel3 = new double[size_c];	

	double * wallDistBFParallel4;
	wallDistBFParallel4 = new double[size_c];	

	double * wallDistAB;
	wallDistAB = new double[size_c];

	double * wallDistAB_t2;
	wallDistAB_t2 = new double[size_c];

	double * wallDistAB_pt1;
	wallDistAB_pt1 = new double[size_c];

	for(int i=0; i<size_c; i++){
		wallDistAB[i] = 1e9;
		wallDistAB_t2[i] = 1e9;
		wallDistAB_pt1[i] = 1e9;
		wallDistBFSerial[i] = 1e9;
		wallDistBFParallel1[i] = 1e9;
		wallDistBFParallel2[i] = 1e9;
		wallDistBFParallel3[i] = 1e9;
		wallDistBFParallel4[i] = 1e9;
	}

	//////////////////////////////////////////////////
	//////////////////////////////////////////////////

	// BRUTEFORCE SERIAL WALLDISTANCE CALCULATION
//	getrusage( RUSAGE_SELF, &tm_start ); 						// Start timer
	SerialBF(size_c,size_f,xc,yc,xf,yf,wallDistBFSerial);
//	getrusage( RUSAGE_SELF, &tm_end );   						// End timer
//	time_in_msec = gettime_msec( &tm_start, &tm_end ); 			// Get elapsed time
	
//	printf("Brute force - serial (main): \t %.0f milliseconds.\n", time_in_msec);
	postproc(ni,nj,wallDistBFSerial,0);							// WRITE TO FILE
	
	// BRUTEFORCE PARALLEL 1 WALLDISTANCE CALCULATION
	getrusage( RUSAGE_SELF, &tm_start ); 						// Start timer
	ParallelBF1(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel1);		
	getrusage( RUSAGE_SELF, &tm_end );   						// End timer
	time_in_msec = gettime_msec( &tm_start, &tm_end ); 			// Get elapsed time

	printf("Brute force - parallel 1 block per cell (main): \t %.0f milliseconds\n\tVerifying output result ...%s\n",
		   time_in_msec, compare_matrices(wallDistBFParallel1, wallDistBFSerial, size_c) ? "Failed" : "Success");

   	// BRUTEFORCE PARALLEL 2 WALLDISTANCE CALCULATION
	getrusage( RUSAGE_SELF, &tm_start ); 						// Start timer
	ParallelBF2(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel2);		
	getrusage( RUSAGE_SELF, &tm_end );   						// End timer
	time_in_msec = gettime_msec( &tm_start, &tm_end ); 			// Get elapsed time

	printf("Brute force - parallel 2 block per cell shared mem (main): \t %.0f milliseconds\n\tVerifying output result ...%s\n",
		   time_in_msec, compare_matrices(wallDistBFParallel2, wallDistBFSerial, size_c) ? "Failed" : "Success");
		   	
	// BRUTEFORCE PARALLEL 3 WALLDISTANCE CALCULATION
	getrusage( RUSAGE_SELF, &tm_start ); 						// Start timer
	ParallelBF3(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel3);		
	getrusage( RUSAGE_SELF, &tm_end );   						// End timer
	time_in_msec = gettime_msec( &tm_start, &tm_end ); 			// Get elapsed time

	printf("Brute force - parallel 3 block per cell coalesced (main): \t %.0f milliseconds\n\tVerifying output result ...%s\n",
		   time_in_msec, compare_matrices(wallDistBFParallel3, wallDistBFSerial, size_c) ? "Failed" : "Success");

   	// BRUTEFORCE PARALLEL 4 WALLDISTANCE CALCULATION
	getrusage( RUSAGE_SELF, &tm_start ); 						// Start timer
	ParallelBF4(size_c,size_f,xc,yc,xf,yf,wallDistBFParallel4);		
	getrusage( RUSAGE_SELF, &tm_end );   						// End timer
	time_in_msec = gettime_msec( &tm_start, &tm_end ); 			// Get elapsed time

	printf("Brute force - parallel 4 thread per cell (main): \t %.0f milliseconds\n\tVerifying output result ...%s\n",
		   time_in_msec, compare_matrices(wallDistBFParallel4, wallDistBFSerial, size_c) ? "Failed" : "Success");


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


	// ABparallel_t1
	getrusage( RUSAGE_SELF, &tm_start ); 						// Start timer
	ab_parallel_t1(xc,yc,xf,yf,size_c,size_f,wallDistAB_pt1);
	getrusage( RUSAGE_SELF, &tm_end );   						// End timer

	postproc(ni,nj,wallDistAB_pt1,3);
	time_in_msec = gettime_msec( &tm_start, &tm_end ); 			// Get elapsed time


	printf("Advancing boundary - parallel T1(main): \t %.0f milliseconds\n", time_in_msec);

	///////////////////////////////////////////////////
	///////////////////////////////////////////////////

	return 0;
}


