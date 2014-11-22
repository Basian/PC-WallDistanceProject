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
#include "timer.h"

#define MAXTHREADSPERBLOCK 512   // This will be used for the thread per cell algorithms

#define MAXTHREADSPERBLOCK1 1024 // Maximum threads per block per devicequery on OSC 
                                 // Probably not good to use, but this may be used in the block per cell 
								 // methods to check algorithm functionality and to benchmark performance


// Kernel for ParallelBF1a: calculates the wall distances for each cell in a block
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

// Brute force - parallel 1a (block per cell method)
void ParallelBF1a(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// 	For timing calculations
	GpuTimer timer;

	// determine the number of threads and blocks
    int threads, blocks;
	
	if (fsize <= MAXTHREADSPERBLOCK1)
	{
	    threads = MAXTHREADSPERBLOCK1;
		blocks = csize;
		while ( (threads/2) >= fsize )
		{
			threads = threads/2;
		}
	}
	else
	{
		// This is currently not supported; in the future, this can be implemented by using a block dimension
		printf("Block per cell method cannot be executed when there are more than %d solid faces\n",MAXTHREADSPERBLOCK1);
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

	timer.Start(); // Start timer
    // Launch block per cell kernel to get the wall distance for each cell, where each cell is a block
    block_per_cell<<<blocks, threads, threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 1a: \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, blocks*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));

}

// Kernel for ParallelBF1b: calculates the wall distances for each cell in a block
__global__ void block_per_cell_shmem(double * dout, int fsize, double * xc, double * yc, double * xf, double * yf)
{
    // sxf, syf, fdistance are allocated in the kernel call: 3rd arg to <<<b, t, shmem>>>
	extern __shared__ double sdata[];
	double *sxf = &sdata[0];			// x of faces
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

// Brute force - parallel 1b (block per cell method using shared memory for optimized reads within a block)
void ParallelBF1b(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// 	For timing calculations
	GpuTimer timer;

	// determine the number of threads and blocks
    int threads, blocks;
	
	if (fsize <= MAXTHREADSPERBLOCK1)
	{
	    threads = MAXTHREADSPERBLOCK1;
		blocks = csize;
		while ( (threads/2) >= fsize )
		{
			threads = threads/2;
		}
	}
	else
	{
		// This is currently not supported; in the future, this can be implemented by using a block dimension
		printf("Block per cell method cannot be executed when there are more than %d solid faces\n",MAXTHREADSPERBLOCK1);
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

	timer.Start(); // Start timer
    // Launch block per cell kernel to get the wall distance for each cell, where each cell is a block
    block_per_cell_shmem<<<blocks, threads, 3 * threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 1b: \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, blocks*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));

}

// Kernel for ParallelBF1c: calculates the wall distances for each cell in a block
__global__ void block_per_cell_shmem_c(double * dout, int fsize, double * xc, double * yc, double * xf, double * yf)
{
    // sxyf, fdistance are allocated in the kernel call: 3rd arg to <<<b, t, shmem>>>
	extern __shared__ double sdata[];
	double *sxyf = &sdata[0]; // This will contain the x and y elements interspersed, i.e. x0, y0, x1, y1, ....
	double *fdistance = &sdata[2*blockDim.x]; 	// Array of distances from each face to a cell center

	// Thread and block IDs
    int tid  = threadIdx.x;
	int bid = blockIdx.x;
	
    // Load shared mem from global mem; sxyf will contain the x and y elements interspersed, i.e. x0, y0, x1, y1, ....
	sxyf[tid*2] = xf[tid];
	sxyf[tid*2 + 1] = yf[tid];	
    __syncthreads();

	// If thread id is less than the number of faces, calculate distance to the cell center
	if ( tid < fsize )
	{
		fdistance[tid] = sqrt(pow(xc[bid]-sxyf[tid*2],2) + pow(yc[bid]-sxyf[tid*2 + 1],2));
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

// Brute force - parallel 1c (block per cell method using shared memory and coalesced reads within a block)
void ParallelBF1c(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// 	For timing calculations
	GpuTimer timer;

	// determine the number of threads and blocks
    int threads, blocks;
	
	if (fsize <= MAXTHREADSPERBLOCK1)
	{
	    threads = MAXTHREADSPERBLOCK1;
		blocks = csize;
		while ( (threads/2) >= fsize )
		{
			threads = threads/2;
		}
	}
	else
	{
		// This is currently not supported; in the future, this can be implemented by using a block dimension
		printf("Block per cell method cannot be executed when there are more than %d solid faces\n",MAXTHREADSPERBLOCK1);
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

	timer.Start(); // Start timer
    // Launch block per cell kernel to get the wall distance for each cell, where each cell is a block
    block_per_cell_shmem_c<<<blocks, threads, 3 * threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 1c: \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, blocks*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));

}


// Kernel for ParallelBF2a: calculates the wall distance for each cell in a thread
__global__ void thread_per_cell(double * dout, int fsize, double * xc, double * yc, double * xf, double * yf)
{
	// Thread and block IDs
	int myid = threadIdx.x + blockDim.x * blockIdx.x;

	double minfdistance = 1e9; // Initialize it to a large value;
	double fdistance;

	// Loop through each face, and calculate its distance from the cell center
	for (int j=0; j<fsize; j++)
	{
//			minfdistance = min(minfdistance,sqrt(pow(xc[myid]-xf[j],2) + pow(yc[myid]-yf[j],2)));
		fdistance = sqrt(pow(xc[myid]-xf[j],2) + pow(yc[myid]-yf[j],2));
		if ( fdistance < minfdistance)
		{
			minfdistance = fdistance;
		}	
	}

    dout[myid] = minfdistance;
	
}

// Brute force - parallel 2a thread per cell method
void ParallelBF2a(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// 	For timing calculations
	GpuTimer timer;

	// determine the number of threads and blocks
    int threads, blocks;
	
    threads = MAXTHREADSPERBLOCK;
	blocks = (csize/threads);
	if(csize%threads) { blocks++; }

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

	timer.Start(); // Start timer
    // Launch thread per cell kernel to get the wall distance for each cell
    thread_per_cell<<<blocks, threads>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 2a: \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, csize*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));

	// DEBUGGING CODE
/*	FILE * fp;
	fp = fopen ("BFParallel4output.txt", "wb\r\n");
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

// Kernel for ParallelBF2b: calculates the wall distance for each cell in a thread (using shared memory)
__global__ void thread_per_cell_sh(double * dout, int fsize, double * xc, double * yc, double * xf, double * yf, int shmem_space)
{
	// Thread and block IDs
	int myid = threadIdx.x + blockDim.x * blockIdx.x;

	double minfdistance = 1e9; // Initialize it to a large value;
	double fdistance;

	// sxf, syf are allocated in the kernel call: 3rd arg to <<<b, t, shmem>>>
	extern __shared__ double sdata[];
	int f_per_i = shmem_space/(2*sizeof(double));
	double *sxf = &sdata[0];		// x of faces
	double *syf = &sdata[f_per_i];	// y of faces
	int loopidx = (2*fsize*sizeof(double))/shmem_space;
	if((2*fsize*sizeof(double))%shmem_space) { loopidx++; }
	
	// Loop through and keep loading allocated shared memory with the correct x and y values
	for (int i=0; i<loopidx; i++)
	{
		// Load shared mem from global mem - do it for one thread only (not good as it causes thread divergence)
		if (threadIdx.x == 0)
		{
			for (int j=0; j<f_per_i; j++){
			if ( i*f_per_i+j < fsize)
				{
				sxf[j] = xf[i*f_per_i+j];
				syf[j] = yf[i*f_per_i+j];
				}
			}
		}
		__syncthreads();

		// Loop through each face, and calculate its distance from the cell center
		for (int j=0; j<f_per_i; j++)
		{
		//	minfdistance = min(minfdistance,sqrt(pow(xc[myid]-xf[j],2) + pow(yc[myid]-yf[j],2)));
			if ( i*f_per_i+j < fsize)
			{
			fdistance = sqrt(pow(xc[myid]-sxf[j],2) + pow(yc[myid]-syf[j],2));
			if ( fdistance < minfdistance)
			{
				minfdistance = fdistance;
			}	
			}
		}
		__syncthreads();
		
	}
    dout[myid] = minfdistance;
	
}

// Brute force - parallel 2b thread per cell method using shared memory
void ParallelBF2b(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// 	For timing calculations
	GpuTimer timer;

	// determine the number of threads and blocks
    int threads, blocks;
	
    threads = MAXTHREADSPERBLOCK;
	blocks = (csize/threads);
	if(csize%threads) { blocks++; }

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
	
	// Calculate shared memory space (max is 49152, but we can limit it to be lower than that)
	int shmem_space = (2 * fsize * sizeof(double) > 49152) ? 49152 : (2 * fsize * sizeof(double));

	timer.Start(); // Start timer
    // Launch thread per cell kernel to get the wall distance for each cell
    thread_per_cell_sh<<<blocks, threads, shmem_space>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf, shmem_space);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 2b: \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, csize*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));


}


// Kernel for ParallelBF2c: calculates the wall distance for each cell in a thread (using coalesced global memory reads)
__global__ void thread_per_cell_c(double * dout, int fsize, double * xc, double * yc, double * xyf)
{
	// Thread and block IDs
	int myid = threadIdx.x + blockDim.x * blockIdx.x;

	double minfdistance = 1e9; // Initialize it to a large value;
	double fdistance;

	// Loop through each face, and calculate its distance from the cell center
	for (int j=0; j<fsize; j++)
	{
//			minfdistance = min(minfdistance,sqrt(pow(xc[myid]-xf[j],2) + pow(yc[myid]-yf[j],2)));
		fdistance = sqrt(pow(xc[myid]-xyf[j*2],2) + pow(yc[myid]-xyf[j*2 + 1],2));
		if ( fdistance < minfdistance)
		{
			minfdistance = fdistance;
		}	
	}

    dout[myid] = minfdistance;
	
}

// Brute force - parallel 2c thread per cell method using shared memory and coalesced reads within a block
void ParallelBF2c(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// 	For timing calculations
	GpuTimer timer;

	// determine the number of threads and blocks
    int threads, blocks;
	
    threads = MAXTHREADSPERBLOCK;
	blocks = (csize/threads);
	if(csize%threads) { blocks++; }

	// declare GPU memory pointers
	double * d_xc, * d_yc, * d_xyf;
    double * d_wallDist;
	double xyf[2*fsize];
	
    // allocate GPU memory
    checkCudaErrors(cudaMalloc(&d_xc, csize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_yc, csize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_xyf, 2*fsize*sizeof(double)));
    checkCudaErrors(cudaMalloc(&d_wallDist, csize*sizeof(double)));

	// Copy from host to GPU memory
    checkCudaErrors(cudaMemcpy(d_xc, xc, csize*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_yc, yc, csize*sizeof(double), cudaMemcpyHostToDevice));
	// For coalesced global memory access, xyf will contain the x and y elements interspersed, i.e. x0, y0, x1, y1, ....
	for ( int j=0; j<fsize; j++)
	{
		xyf[j*2] = xf[j];
		xyf[j*2+1] = yf[j];
	}
	checkCudaErrors(cudaMemcpy(d_xyf, xyf, 2*fsize*sizeof(double), cudaMemcpyHostToDevice));

	timer.Start(); // Start timer
    // Launch thread per cell kernel to get the wall distance for each cell
 	thread_per_cell_c<<<blocks, threads>>>(d_wallDist, fsize, d_xc, d_yc, d_xyf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 2c: \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, csize*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xyf));
    checkCudaErrors(cudaFree(d_wallDist));


}
