/*
 * ABserial.c
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#include "ab_pt1.h"

extern "C" {
#include "boundBox.h"
}

extern "C" {
#include "computeAuxiliaryGrid_pt1.h"
}

extern "C" {
#include "compactAuxiliaryGrid_pt1.h"
}
//#include "writecell.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <sys/time.h>
#include "cuda_utils.h"
#include "timer.h"
#include <time.h>

__global__
void wd_ab_parallel_t1(double *xc, double *yc, double *xf, double *yf, double *xbox, double *ybox, struct cell_pt1 *compAuxCells, int size_c, int size_f, int numAuxCells, double auxDiag, double *wallDist){


	int myId, includesAuxCells, j, index;
	double r, rtemp, rAux;

	myId = threadIdx.x + blockDim.x * blockIdx.x;

	// Keep array access bounded
	if (myId >= size_c){
		return;
	}

	// Compute initial radius
	r=1e9;
	for (j=0; j<8; j++){
		rtemp = sqrt( pow((xc[myId]-xbox[j]),2) + pow((yc[myId]-ybox[j]),2) );
		if (rtemp<r){
			r=rtemp;
		}
	}



	// Loop through compacted auxCell array to see if any lie within rc
	includesAuxCells = 0;
	while(includesAuxCells == 0){
		for (j=0; j<numAuxCells; j++){

			rAux = sqrt( pow(xc[myId]-compAuxCells[j].xcenter,2) + pow(yc[myId]-compAuxCells[j].ycenter,2) );
			// Increase r to be sure enough geometry is included
			if(rAux < r){
				r += auxDiag*0.5;
				includesAuxCells=1;
				break;
			}
			else{
				r += auxDiag;
			}
//			if(myId==0){
//				printf("rAux, r: %f, %f\n",rAux,r);
//			}
		}

	}


	/*
	 *  Loop through compacted auxCell array. For those that lie within r,
	 *  traverse through faces, compute wallDist and check for minimum
	 */
	for (j=0; j<numAuxCells; j++){

		rAux = sqrt( pow(xc[myId]-compAuxCells[j].xcenter,2) + pow(yc[myId]-compAuxCells[j].ycenter,2));

		// Check if auxCell is within radius of interest
		if(rAux < r){
			index = 0;

			// Loop through faces and compute distance from grid cell center
			while(index < compAuxCells[j].numFaces){
				rtemp = sqrt( pow(xc[myId]-xf[compAuxCells[j].faceIndex[index]],2) + pow(yc[myId]-yf[compAuxCells[j].faceIndex[index]],2));

				// If dist is smaller than current wallDist, replace
				if(rtemp<wallDist[myId]){
					wallDist[myId]=rtemp;
				}

				index++;
			}
		}

	}





}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


__global__
void wd_ab_parallel_t2(double *cellCenters, double *xf, double *yf, double *xbox, double *ybox, struct cell_pt1 *compAuxCells, int size_c, int size_f, int numAuxCells, double auxDiag, double *wallDist){

	extern __shared__ double s_cellCenters[];


	int myId, tid, includesAuxCells, j, index;
	double r, rtemp, rAux;

	tid = threadIdx.x;
	myId = threadIdx.x + blockDim.x * blockIdx.x;


	// Pull cell centers into shared memory
	s_cellCenters[2*tid] = cellCenters[2*myId];
	s_cellCenters[2*tid+1] = cellCenters[2*myId+1];



	// Keep array access bounded
	if (myId >= size_c){
		return;
	}

	// Compute initial radius
	r=1e9;
	for (j=0; j<8; j++){
		rtemp = sqrt( pow((s_cellCenters[2*tid]-xbox[j]),2) + pow((s_cellCenters[2*tid+1]-ybox[j]),2) );
		if (rtemp<r){
			r=rtemp;
		}
	}



	// Loop through compacted auxCell array to see if any lie within rc
	includesAuxCells = 0;
	while(includesAuxCells == 0){
		for (j=0; j<numAuxCells; j++){

			rAux = sqrt( pow(s_cellCenters[2*tid]-compAuxCells[j].xcenter,2) + pow(s_cellCenters[2*tid+1]-compAuxCells[j].ycenter,2) );
			// Increase r to be sure enough geometry is included
			if(rAux < r){
				r += auxDiag*0.5;
				includesAuxCells=1;
				break;
			}
			else{
				r += auxDiag;
			}
//			if(myId==0){
//				printf("rAux, r: %f, %f\n",rAux,r);
//			}
		}

	}


	/*
	 *  Loop through compacted auxCell array. For those that lie within r,
	 *  traverse through faces, compute wallDist and check for minimum
	 */
	for (j=0; j<numAuxCells; j++){

		rAux = sqrt( pow(s_cellCenters[2*tid]-compAuxCells[j].xcenter,2) + pow(s_cellCenters[2*tid+1]-compAuxCells[j].ycenter,2));

		// Check if auxCell is within radius of interest
		if(rAux < r){
			index = 0;

			// Loop through faces and compute distance from grid cell center
			while(index < compAuxCells[j].numFaces){
				rtemp = sqrt( pow(s_cellCenters[2*tid]-xf[compAuxCells[j].faceIndex[index]],2) + pow(s_cellCenters[2*tid+1]-yf[compAuxCells[j].faceIndex[index]],2));

				// If dist is smaller than current wallDist, replace
				if(rtemp<wallDist[myId]){
					wallDist[myId]=rtemp;
				}

				index++;
			}
		}

	}





}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



__global__
void wd_ab_parallel_t3(double *cellCenters, double *faceCenters, double *xbox, double *ybox, struct cell_pt1 *compAuxCells, int size_c, int size_f, int numAuxCells, double auxDiag, double *wallDist){

	extern __shared__ double s_faceCenters[];


	int myId, tid, includesAuxCells, i, j, index;
	double r, rtemp, rcurrent, rAux;

	tid = threadIdx.x;
	myId = threadIdx.x + blockDim.x * blockIdx.x;


	// Pull face centers into shared memory
	if (tid==0){
		for (i=0; i<size_f; i++){
			s_faceCenters[2*i] = faceCenters[2*i];
			s_faceCenters[2*i+1] = faceCenters[2*i+1];
		}
	}
	__syncthreads();

	// Keep array access bounded
	if (myId >= size_c){
		return;
	}

	// Compute initial radius
	r=1e9;
	for (j=0; j<8; j++){
		rtemp = sqrt( pow((cellCenters[2*myId]-xbox[j]),2) + pow((cellCenters[2*myId+1]-ybox[j]),2) );
		if (rtemp<r){
			r=rtemp;
		}
	}



	// Loop through compacted auxCell array to see if any lie within rc
	includesAuxCells = 0;
	while(includesAuxCells == 0){
		for (j=0; j<numAuxCells; j++){

			rAux = sqrt( pow(cellCenters[2*myId]-compAuxCells[j].xcenter,2) + pow(cellCenters[2*myId+1]-compAuxCells[j].ycenter,2) );
			// Increase r to be sure enough geometry is included
			if(rAux < r){
				r += auxDiag*0.5;
				includesAuxCells=1;
				break;
			}
			else{
				r += auxDiag;
			}
//			if(myId==0){
//				printf("rAux, r: %f, %f\n",rAux,r);
//			}
		}

	}


	/*
	 *  Loop through compacted auxCell array. For those that lie within r,
	 *  traverse through faces, compute wallDist and check for minimum
	 */
	for (j=0; j<numAuxCells; j++){

		rAux = sqrt( pow(cellCenters[2*myId]-compAuxCells[j].xcenter,2) + pow(cellCenters[2*myId+1]-compAuxCells[j].ycenter,2));

		// Check if auxCell is within radius of interest
		if(rAux < r){
			index = 0;

			// Loop through faces and compute distance from grid cell center
			while(index < compAuxCells[j].numFaces){
				rtemp = sqrt( pow(cellCenters[2*myId]-s_faceCenters[compAuxCells[j].faceIndex[index]],2) + pow(cellCenters[2*myId+1]-s_faceCenters[compAuxCells[j].faceIndex[index]],2));

				// If dist is smaller than current wallDist, replace
				if(rtemp<rcurrent){
//					wallDist[myId]=rtemp;
					rcurrent = rtemp;

				}

				index++;
			}
		}

	}

	// Store wallDistance to global array
	wallDist[myId] = rcurrent;


}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


__global__
void wd_ab_parallel_t4(double *cellCenters, double *faceCenters, double *xbox, double *ybox, struct cell_pt1 *compAuxCells, int size_c, int size_f, int numAuxCells, double auxDiag, double *wallDist){



	int myId, tid, includesAuxCells, i, j, index;
	double r, rtemp, rcurrent, rAux, c_x, c_y;

	tid = threadIdx.x;
	myId = threadIdx.x + blockDim.x * blockIdx.x;


	// Keep array access bounded
	if (myId >= size_c){
		return;
	}


	c_x = cellCenters[2*myId];
	c_y = cellCenters[2*myId+1];



	// Compute initial radius
	r=1e9;
	for (j=0; j<8; j++){
		rtemp = sqrt( pow((c_x-xbox[j]),2) + pow((c_y-ybox[j]),2) );
		if (rtemp<r){
			r=rtemp;
		}
	}



	// Loop through compacted auxCell array to see if any lie within rc
	includesAuxCells = 0;
	while(includesAuxCells == 0){
		for (j=0; j<numAuxCells; j++){

			rAux = sqrt( pow(c_x-compAuxCells[j].xcenter,2) + pow(c_y-compAuxCells[j].ycenter,2) );
			// Increase r to be sure enough geometry is included
			if(rAux < r){
				r += auxDiag*0.5;
				includesAuxCells=1;
				break;
			}
			else{
				r += auxDiag;
			}

		}

	}


	/*
	 *  Loop through compacted auxCell array. For those that lie within r,
	 *  traverse through faces, compute wallDist and check for minimum
	 */
	for (j=0; j<numAuxCells; j++){

		rAux = sqrt( pow(c_x-compAuxCells[j].xcenter,2) + pow(c_y-compAuxCells[j].ycenter,2));

		// Check if auxCell is within radius of interest
		if(rAux < r){
			index = 0;

			// Loop through faces and compute distance from grid cell center
			while(index < compAuxCells[j].numFaces){
				rtemp = sqrt( pow(c_x-faceCenters[compAuxCells[j].faceIndex[index]],2) + pow(c_y-faceCenters[compAuxCells[j].faceIndex[index]],2));

				// If dist is smaller than current wallDist, replace
				if(rtemp<rcurrent){
//					wallDist[myId]=rtemp;
					rcurrent = rtemp;

				}

				index++;
			}
		}

	}

	// Store wallDistance to global array
	wallDist[myId] = rcurrent;


}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


__global__
void wd_ab_parallel_t5(double *cellCenters, double *faceCenters, double *xbox, double *ybox, struct cell_pt1 *compAuxCells, int size_c, int size_f, int numAuxCells, double auxDiag, double *wallDist){

	extern __shared__ double s_faceCenters [];


	int myId, tid, includesAuxCells, i, j, index, f_start, f_end, s_ind;
	double r, rtemp, rcurrent, rAux, c_x, c_y;

	tid = threadIdx.x;
	myId = threadIdx.x + blockDim.x * blockIdx.x;


	// Keep array access bounded
	if (myId >= size_c){
		return;
	}

	f_start=0;
	f_end=1024;


	// Loop through shared memory cycles until all faces have been processed
//	while


		// Pull 1024 faces into shared memory
		if (tid==0){
			s_ind = 0;
			for (i=f_start; i<f_end; i++){
				s_faceCenters[2*i] = faceCenters[2*myId];
				s_faceCenters[2*i+1] = faceCenters[2*myId+1];
			}
		}

	c_x = cellCenters[2*myId];
	c_y = cellCenters[2*myId+1];



	// Compute initial radius
	r=1e9;
	for (j=0; j<8; j++){
		rtemp = sqrt( pow((c_x-xbox[j]),2) + pow((c_y-ybox[j]),2) );
		if (rtemp<r){
			r=rtemp;
		}
	}



	// Loop through compacted auxCell array to see if any lie within rc
	includesAuxCells = 0;
	while(includesAuxCells == 0){
		for (j=0; j<numAuxCells; j++){

			rAux = sqrt( pow(c_x-compAuxCells[j].xcenter,2) + pow(c_y-compAuxCells[j].ycenter,2) );
			// Increase r to be sure enough geometry is included
			if(rAux < r){
				r += auxDiag*0.5;
				includesAuxCells=1;
				break;
			}
			else{
				r += auxDiag;
			}

		}

	}


	/*
	 *  Loop through compacted auxCell array. For those that lie within r,
	 *  traverse through faces, compute wallDist and check for minimum
	 */
	for (j=0; j<numAuxCells; j++){

		rAux = sqrt( pow(c_x-compAuxCells[j].xcenter,2) + pow(c_y-compAuxCells[j].ycenter,2));

		// Check if auxCell is within radius of interest
		if(rAux < r){
			index = 0;

			// Loop through faces and compute distance from grid cell center
			while(index < compAuxCells[j].numFaces){
				rtemp = sqrt( pow(c_x-faceCenters[compAuxCells[j].faceIndex[index]],2) + pow(c_y-faceCenters[compAuxCells[j].faceIndex[index]],2));

				// If dist is smaller than current wallDist, replace
				if(rtemp<rcurrent){
//					wallDist[myId]=rtemp;
					rcurrent = rtemp;

				}

				index++;
			}
		}

	}

	// Store wallDistance to global array
	wallDist[myId] = rcurrent;


}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



void ab_parallel_t1(double * xc, double * yc, double * xf, double * yf, int size_c, int size_f, double * wallDist){

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
	int resI=80;
	int resJ=80;
	double auxDiag = sqrt( pow((xmax-xmin)/(double)(resI-1),2) + pow((ymax-ymin)/(double)(resJ-1),2));


	int numAuxCells = (resI-1)*(resJ-1);
	int i, j, cellsWithFaces;
	struct cell_pt1 *auxCells;
//	auxCells = (struct cell_pt1 *)malloc(numAuxCells*sizeof(struct cell_pt1));
	auxCells = new cell_pt1[numAuxCells];


	computeAuxiliaryGrid_pt1(xmin,xmax,ymin,ymax,resI,resJ,auxCells);

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
	struct cell_pt1 * compAuxCells;
//	compAuxCells = (struct cell_pt1 *)malloc(cellsWithFaces*sizeof(struct cell_pt1));
	compAuxCells = new cell_pt1[cellsWithFaces];


	///////
	compactAuxiliaryGrid_pt1(auxCells,numAuxCells,compAuxCells,xf,yf,size_f);
	///////




	// Bounding box point arrays
	double xmid = (xmax+xmin)/2.0;
	double ymid = (ymax+ymin)/2.0;

	double xBoxPts[8] = {xmin, xmid, xmax, xmax, xmax, xmid, xmin, xmin};
	double yBoxPts[8] = {ymin, ymin, ymin, ymid, ymax, ymax, ymax, ymid};



	////////////////////////////////////////////////////////////////////////////////
	//	Combine xc,yc arrays for coallesced memory access in parallel t2 version
	////////////////////////////////////////////////////////////////////////////////

	double *cellCenters;
	cellCenters = new double[2*size_c];

	for (i=0; i<size_c; i++){
		cellCenters[2*i] = xc[i];
		cellCenters[2*i+1] = yc[i];
	}


	double *faceCenters;
	faceCenters = new double[2*size_f];

	for (i=0; i<size_f; i++){
		faceCenters[2*i] = xf[i];
		faceCenters[2*i+1] = yf[i];
	}



	////////////////////////////////////////////////////////////////////
	//  Allocate device memory and copy data
	////////////////////////////////////////////////////////////////////
	// bounding box
	double *d_xbox, *d_ybox;
	checkCudaErrors(cudaMalloc(&d_xbox,8*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_ybox,8*sizeof(double)));

	checkCudaErrors(cudaMemcpy(d_xbox,xBoxPts,8*sizeof(double),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_ybox,yBoxPts,8*sizeof(double),cudaMemcpyHostToDevice));

	// grid and faces
	double *d_xc, *d_yc, *d_xf, *d_yf, *d_cellCenters, *d_faceCenters;
	checkCudaErrors(cudaMalloc(&d_xc,size_c*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_yc,size_c*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_xf,size_c*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_yf,size_c*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_cellCenters,2*size_c*sizeof(double)));
	checkCudaErrors(cudaMalloc(&d_faceCenters,2*size_f*sizeof(double)));

	checkCudaErrors(cudaMemcpy(d_xc,xc,size_c*sizeof(double),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_yc,yc,size_c*sizeof(double),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_xf,xf,size_c*sizeof(double),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_yf,yf,size_c*sizeof(double),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_cellCenters,cellCenters,2*size_c*sizeof(double),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_faceCenters,faceCenters,2*size_f*sizeof(double),cudaMemcpyHostToDevice));

	// auxCell structs
	struct cell_pt1 * d_compAuxCells;
	checkCudaErrors(cudaMalloc((void **)&d_compAuxCells,cellsWithFaces*sizeof(struct cell_pt1)));
	checkCudaErrors(cudaMemcpy(d_compAuxCells,compAuxCells,cellsWithFaces*sizeof(struct cell_pt1),cudaMemcpyHostToDevice));


	// wallDist array
	double *d_wallDist;
	checkCudaErrors(cudaMalloc(&d_wallDist,size_c*sizeof(double)));
	checkCudaErrors(cudaMemcpy(d_wallDist,wallDist,size_c*sizeof(double),cudaMemcpyHostToDevice));




	////////////////////////////////////////////////////////////////////
	//	Wall Distance Calc
	////////////////////////////////////////////////////////////////////
	GpuTimer timer;

	int threadsPerBlock, numBlocks;
	threadsPerBlock = 512;
	numBlocks = (size_c/threadsPerBlock)+1;



	timer.Start();
	wd_ab_parallel_t1<<<numBlocks,threadsPerBlock>>>(d_xc,d_yc,d_xf,d_yf,d_xbox,d_ybox,d_compAuxCells,size_c,size_f,cellsWithFaces,auxDiag,d_wallDist);
	timer.Stop();
	printf("Advancing boundary - parallel T1(GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());


	// Reset wallDistance
	checkCudaErrors(cudaMemcpy(d_wallDist,wallDist,size_c*sizeof(double),cudaMemcpyHostToDevice));



	timer.Start();
	wd_ab_parallel_t2<<<numBlocks,threadsPerBlock,2*threadsPerBlock*sizeof(double)>>>(d_cellCenters,d_xf,d_yf,d_xbox,d_ybox,d_compAuxCells,size_c,size_f,cellsWithFaces,auxDiag,d_wallDist);
	timer.Stop();
	printf("Advancing boundary - parallel T2(GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());


//	timer.Start();
//	wd_ab_parallel_t3<<<numBlocks,threadsPerBlock,2*size_f*sizeof(double)>>>(d_cellCenters,d_faceCenters,d_xbox,d_ybox,d_compAuxCells,size_c,size_f,cellsWithFaces,auxDiag,d_wallDist);
//	timer.Stop();
//	printf("Advancing boundary - parallel T3(GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());


	timer.Start();
	wd_ab_parallel_t4<<<numBlocks,threadsPerBlock>>>(d_cellCenters,d_faceCenters,d_xbox,d_ybox,d_compAuxCells,size_c,size_f,cellsWithFaces,auxDiag,d_wallDist);
	timer.Stop();
	printf("Advancing boundary - parallel T4(GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());


//	timer.Start();
//	wd_ab_parallel_t5<<<numBlocks,threadsPerBlock,2048*sizeof(double)>>>(d_cellCenters,d_faceCenters,d_xbox,d_ybox,d_compAuxCells,size_c,size_f,cellsWithFaces,auxDiag,d_wallDist);
//	timer.Stop();
//	printf("Advancing boundary - parallel T5(GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());



	// Copy wallDist back to host
	checkCudaErrors(cudaMemcpy(wallDist,d_wallDist,sizeof(double)*size_c,cudaMemcpyDeviceToHost));







	////////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////////



}




