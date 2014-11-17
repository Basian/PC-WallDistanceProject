/*
 * ABserial.c
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#include "ab_serial_t2.h"
#include "boundBox.h"
#include "computeAuxiliaryGrid_t2.h"
#include "compactAuxiliaryGrid_t2.h"
//#include "writecell.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <sys/times.h>

void ab_serial_t2(double * xc, double * yc, double * xf, double * yf, int size_c, int size_f, double * wallDist){

	double xmin;
	double xmax;
	double ymin;
	double ymax;

	////////////////////////////////////////////////////////////////////
	//		Pre-processing
	////////////////////////////////////////////////////////////////////

	// Create geometry bounding box
	boundBox(xf,yf,size_f,&xmin,&xmax,&ymin,&ymax);


	// Create auxiliary grid
	int resI=20;
	int resJ=20;
	double auxDiag = sqrt( pow((xmax-xmin)/(double)(resI-1),2) + pow((ymax-ymin)/(double)(resJ-1),2));
	int numAuxCells = (resI-1)*(resJ-1);
	int i, j, k, numFaces, cellsWithFaces;
	struct cell_t2 * auxCells;
	auxCells = (struct cell_t2 *)malloc(numAuxCells*sizeof(struct cell_t2));

	computeAuxiliaryGrid_t2(xmin,xmax,ymin,ymax,resI,resJ,auxCells);
//	writecell(auxCells,numAuxCells,0);

	// Count number of auxiliary cells that contain geometry faces
	cellsWithFaces = 0;
	for (i=0; i<numAuxCells; i++){
		for (j=0; j<size_f; j++){

			if (xf[j] < auxCells[i].xmax && xf[j] > auxCells[i].xmin && yf[j] < auxCells[i].ymax && yf[j] > auxCells[i].ymin){
				cellsWithFaces++;
				break;
			}

		}
	}

	// Allocate memory for compacted cells
	struct cell_t2 * compAuxCells;
	compAuxCells = (struct cell_t2 *)malloc(cellsWithFaces*sizeof(struct cell_t2));

	// Initialize face arrays
	for (i=0; i<numAuxCells; i++){
		auxCells[i].storage = 2;
		auxCells[i].xface = (double *)malloc(2*sizeof(double));
		auxCells[i].yface = (double *)malloc(2*sizeof(double));
	}


	compactAuxiliaryGrid_t2(auxCells,numAuxCells,compAuxCells,xf,yf,size_f);

//	writecell(compAuxCells,cellsWithFaces,1);
//	writefaces(compAuxCells,cellsWithFaces);

	////////////////////////////////////////////////////////////////////
	//	Wall Distance Calc
	////////////////////////////////////////////////////////////////////
	struct rusage tm_start;
	struct rusage tm_end;
	getrusage( RUSAGE_SELF, &tm_start );


	/*
	 * For each grid cell, set an initial radius and gradually increase until
	 * come auxiliary cells lie within the radius. Search the faces included
	 * in those auxiliary cells
	 */
	int index;
	double rc, rAux, rFace;
	double xmid = (xmax+xmin)/2.0;
	double ymid = (ymax+ymin)/2.0;

	double xBoxPts[8] = {xmin, xmid, xmax, xmax, xmax, xmid, xmin, xmin};
	double yBoxPts[8] = {ymin, ymin, ymin, ymid, ymax, ymax, ymax, ymid};
	double rBoxPts[8];

	// Loop through grid cells
	for (i=0; i<size_c; i++){

		// Compute radius of 8 box points
		for (j=0; j<8; j++){
			rBoxPts[j] = sqrt( pow((xc[i]-xBoxPts[j]),2) + pow((yc[i]-yBoxPts[j]),2) );
		}

		// Find minimum of box point radii
		rc = rBoxPts[0];
		for (j=1; j<8; j++){
			if (rBoxPts[j] < rc){
				rc = rBoxPts[j];
			}
		}

		// Loop through compacted auxCell array to see if any lie within rc
		int iterationsA, iterationsB;
		iterationsA=0;
		iterationsB=0;

		int includesAuxCells = 0;
		while(includesAuxCells == 0){
			for (j=0; j<cellsWithFaces; j++){

				rAux = sqrt( pow(xc[i]-compAuxCells[j].xcenter,2) + pow(yc[i]-compAuxCells[j].ycenter,2) );
				// Increase rc to be sure enough geometry is included
				if(rAux < rc){
					rc += auxDiag*0.5;
					includesAuxCells=1;
					break;
				}
				else{
					rc += auxDiag;
				}
//				iterationsA++;

			}
//			iterationsB++;
		}
//		printf("Iterations (A/B): %i, %i\n",iterationsA,iterationsB);

		/*
		 *  Loop through compacted auxCell array. For those that lie within rc,
		 *  traverse through faces, compute wallDist and check for minimum
		 */
		for (j=0; j<cellsWithFaces; j++){

			rAux = sqrt( pow(xc[i]-compAuxCells[j].xcenter,2) + pow(yc[i]-compAuxCells[j].ycenter,2));

			// Check if auxCell is within radius of interest
			if(rAux < rc){
				index = 0;

				// Loop through faces
				while(index < compAuxCells[j].faceNum){
					rFace = sqrt( pow(xc[i]-compAuxCells[j].xface[index],2) + pow(yc[i]-compAuxCells[j].yface[index],2));
					if(rFace<wallDist[i]){
						wallDist[i]=rFace;
					}
					index++;
				}
			}

		}


	}


	getrusage( RUSAGE_SELF, &tm_end );

	double end = (double)tm_end.ru_utime.tv_sec + (double)tm_end.ru_utime.tv_usec / 1000000.0;
	double start = (double)tm_start.ru_utime.tv_sec + (double)tm_start.ru_utime.tv_usec / 1000000.0;

	double diff = end-start;
	printf("Advancing boundary - serial T2: \t \t %.0f milliseconds\n", diff*1000);


	////////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////////



}




