/*
 * bruteforce_parallel.cu
 *
 *  Created on: Nov 15, 2014
 *      Author: vasanth
 */

#include <math.h>
#include <stdio.h>
#include "bruteforce_parallel.h"
#include <time.h>
#include <cuda.h>
#include "cuda_utils.h"
#include <sys/resource.h>
#include <sys/times.h>

#define MAXTHREADSPERBLOCK 1024

// Kernel for ParallelBF1: calculates the wall distances for each cell in a block
__global__ void block_per_cell(double * dout, int fsize, double * xc, double * yc, double * xf, double * yf)
{
    // fdistance is allocated in the kernel call: 3rd arg to <<<b, t, shmem>>>
	// Array of distances from each face to a cell center
    extern __shared__ double fdistance[];
	
	// Thread and block IDs
    int tid  = threadIdx.x;
	int bid = blockIdx.x;

	// If thread id is less than the number of faces, calculate distance to the cell center
	if ( tid < fsize )
	{
		fdistance[tid] = sqrt(pow(xc[bid]-xf[tid],2) + pow(yc[bid]-yf[tid],2));
	}
	// Else, initialize it to a large value; this is necessary for the reduce step to work correctly
	else
	{
		fdistance[tid] = 1e9;
	}
	
	__syncthreads();            // make sure that all the distances in the entire block are calculated!
	
    // Reduce min
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            fdistance[tid] = min(fdistance[tid], fdistance[tid + s]);
        }
        __syncthreads();        // make sure all are done before going to next iteration!
    }

    // Only thread 0 writes result for this block back to global memory
    if (tid == 0)
    {
        dout[bid] = fdistance[0];
    }
}

// Brute force - parallel 1 (block per cell method)
void ParallelBF1(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// determine the number of threads and blocks
    int threads, blocks;
	
	if (fsize <= MAXTHREADSPERBLOCK)
	{
	    threads = MAXTHREADSPERBLOCK;
		blocks = csize;
		while ( (threads/2) >= fsize )
		{
			threads = threads/2;
		}
	}
	else
	{
		printf("Block per cell method cannot be executed when there are more than %d solid faces",MAXTHREADSPERBLOCK);
		return;
	}

	// declare GPU memory pointers
	double * d_xc, * d_yc, * d_xf, * d_yf;
    double * d_wallDist;
	
    // allocate GPU memory
    checkCudaErrors(cudaMalloc(&d_xc, csize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_yc, csize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_xf, fsize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_yf, fsize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_wallDist, csize*sizeof(double)));

	// Copy from host to GPU memory
    checkCudaErrors(cudaMemcpy(d_xc, xc, csize*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_yc, yc, csize*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_xf, xf, fsize*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_yf, yf, fsize*sizeof(double), cudaMemcpyHostToDevice));

    // Launch block per cell kernel to get the wall distance for each cell, where each cell is a block
    block_per_cell<<<blocks, threads, threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, blocks*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));

}

// Kernel for ParallelBF2: calculates the wall distances for each cell in a block
__global__ void block_per_cell_shmem(double * dout, int fsize, double * xc, double * yc, double * xf, double * yf)
{
    // sxf, syf, fdistance are allocated in the kernel call: 3rd arg to <<<b, t, shmem>>>
	extern __shared__ double sdata[];
	double *sxf = &sdata[0];;			// x of faces
	double *syf = &sdata[blockDim.x];	// y of faces
	double *fdistance = &sdata[2*blockDim.x]; 	// Array of distances from each face to a cell center

	// Thread and block IDs
    int tid  = threadIdx.x;
	int bid = blockIdx.x;
	
    // Load shared mem from global mem
    sxf[tid] = xf[tid];
	syf[tid] = yf[tid];
    __syncthreads();

	// If thread id is less than the number of faces, calculate distance to the cell center
	if ( tid < fsize )
	{
		fdistance[tid] = sqrt(pow(xc[bid]-sxf[tid],2) + pow(yc[bid]-syf[tid],2));
	}
	// Else, initialize it to a large value; this is necessary for the reduce step to work correctly
	else
	{
		fdistance[tid] = 1e9;
	}
	
	__syncthreads();            // make sure that all the distances in the entire block are calculated!
	
    // Reduce min
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            fdistance[tid] = min(fdistance[tid], fdistance[tid + s]);
        }
        __syncthreads();        // make sure all are done before going to next iteration!
    }

    // Only thread 0 writes result for this block back to global memory
    if (tid == 0)
    {
        dout[bid] = fdistance[0];
    }
}

// Brute force - parallel 2 (block per cell method using shared memory for optimized reads within a block)
void ParallelBF2(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// determine the number of threads and blocks
    int threads, blocks;
	
	if (fsize <= MAXTHREADSPERBLOCK)
	{
	    threads = MAXTHREADSPERBLOCK;
		blocks = csize;
		while ( (threads/2) >= fsize )
		{
			threads = threads/2;
		}
	}
	else
	{
		printf("Block per cell method cannot be executed when there are more than %d solid faces",MAXTHREADSPERBLOCK);
		return;
	}
	
	// declare GPU memory pointers
	double * d_xc, * d_yc, * d_xf, * d_yf;
    double * d_wallDist;
	
    // allocate GPU memory
    checkCudaErrors(cudaMalloc(&d_xc, csize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_yc, csize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_xf, fsize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_yf, fsize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_wallDist, csize*sizeof(double)));

	// Copy from host to GPU memory
    checkCudaErrors(cudaMemcpy(d_xc, xc, csize*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_yc, yc, csize*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_xf, xf, fsize*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_yf, yf, fsize*sizeof(double), cudaMemcpyHostToDevice));

    // Launch block per cell kernel to get the wall distance for each cell, where each cell is a block
    block_per_cell_shmem<<<blocks, threads, 3 * threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, blocks*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));

}
