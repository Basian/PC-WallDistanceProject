/*
 * bruteforce.cxx
 *
 *  Created on: Nov 11, 2014
 *      Author: vasanth
 */

#include <math.h>
#include <stdio.h>
#include "bruteforce.h"
#include <sys/resource.h>
#include <sys/times.h>

void SerialBF(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist){

	// Used for timing calculations
	struct rusage tm_start;
	struct rusage tm_end;

	// Array of distances from each face to a cell center
	double fdistance[fsize];
	// Minimum face distance ("running" minimum)
	double minfdistance;

	getrusage( RUSAGE_SELF, &tm_start ); // Start timer
	
	// Loop through each cell center and find its wall distance as the minimum distance from the solid faces
	for (int i=0; i<csize; i++){
		// Loop through each face, and calculate its distance from the cell center
		for (int j=0; j<fsize; j++)
		{
			fdistance[j] = sqrt(pow(xc[i]-xf[j],2) + pow(yc[i]-yf[j],2));
			// Initialize minimum face distance in first iteration
			if ( j == 0 )
			{
				minfdistance = fdistance[j];
			}
			// Else, update minimum face distance if the face distance is smaller than the current minimum
			else if ( fdistance[j] < minfdistance)
			{
				minfdistance = fdistance[j];
			}			
		}
		// The cell's wall distance is the minimum distance from the solid faces
		wallDist[i] = minfdistance;
	}
	
	// End timer and calculate time
	getrusage( RUSAGE_SELF, &tm_end );

	double end = (double)tm_end.ru_utime.tv_sec + (double)tm_end.ru_utime.tv_usec / 1000000.0;
	double start = (double)tm_start.ru_utime.tv_sec + (double)tm_start.ru_utime.tv_usec / 1000000.0;

	double diff = end-start;
	printf("Brute force - serial: \t \t %.0f milliseconds\n", diff*1000);

/* 	// DEBUGGING CODE
	FILE * fp;
	fp = fopen ("BFSerialoutput.txt", "wb\r\n");
	fprintf(fp, "Cell center (x), Cell center (y), Wall Distance \n");
	for (int i = 0; i < csize; i++ )
	{
		fprintf(fp, "%lf, %lf, %lf\n", xc[i], yc[i], wallDist[i]);
	}
	fprintf(fp, "Face center (x), Face center (y)\n");
	for (int i = 0; i < fsize; i++ )
	{
		fprintf(fp, "%lf, %lf\n", xf[i], yf[i]);
	}  
	fclose(fp);
 */
}
