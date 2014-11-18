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
#include "timer.h"

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
	// 	For timing calculations
	GpuTimer timer;

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
		// This is currently not supported; in the future, this can be implemented by using a block dimension
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

	timer.Start(); // Start timer
    // Launch block per cell kernel to get the wall distance for each cell, where each cell is a block
    block_per_cell<<<blocks, threads, threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 1 block per cell (GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());

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

// Brute force - parallel 2 (block per cell method using shared memory for optimized reads within a block)
void ParallelBF2(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// 	For timing calculations
	GpuTimer timer;

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
		// This is currently not supported; in the future, this can be implemented by using a block dimension
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

	timer.Start(); // Start timer
    // Launch block per cell kernel to get the wall distance for each cell, where each cell is a block
    block_per_cell_shmem<<<blocks, threads, 3 * threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 2 block per cell shared mem (GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, blocks*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));

}

// Kernel for ParallelBF3: calculates the wall distances for each cell in a block
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

// Brute force - parallel 3 (block per cell method using shared memory and coalesced reads within a block)
void ParallelBF3(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
{
	// 	For timing calculations
	GpuTimer timer;

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
		// This is currently not supported; in the future, this can be implemented by using a block dimension
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

	timer.Start(); // Start timer
    // Launch block per cell kernel to get the wall distance for each cell, where each cell is a block
    block_per_cell_shmem_c<<<blocks, threads, 3 * threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 3 block per cell coalesced (GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, blocks*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));

}


// Kernel for ParallelBF4: calculates the wall distance for each cell in a thread
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

// Brute force - parallel 4 thread per cell method
void ParallelBF4(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
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
	printf("Brute force - parallel 4 thread per cell (GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());

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

// Kernel for ParallelBF5: calculates the wall distance for each cell in a thread (using shared memory)
__global__ void thread_per_cell_sh(double * dout, int fsize, double * xc, double * yc, double * xf, double * yf)
{
	// Thread and block IDs
	int myid = threadIdx.x + blockDim.x * blockIdx.x;
    int tid  = threadIdx.x;

	// sxf, syf are allocated in the kernel call: 3rd arg to <<<b, t, shmem>>>
	extern __shared__ double sdata[];
	double *sxf = &sdata[0];			// x of faces
	double *syf = &sdata[blockDim.x];	// y of faces
	
    // Load shared mem from global mem
    sxf[tid] = xf[tid];
	syf[tid] = yf[tid];
    __syncthreads();

	double minfdistance = 1e9; // Initialize it to a large value;
	double fdistance;

	// Loop through each face, and calculate its distance from the cell center
	for (int j=0; j<fsize; j++)
	{
//			minfdistance = min(minfdistance,sqrt(pow(xc[myid]-xf[j],2) + pow(yc[myid]-yf[j],2)));
		fdistance = sqrt(pow(xc[myid]-sxf[j],2) + pow(yc[myid]-syf[j],2));
		if ( fdistance < minfdistance)
		{
			minfdistance = fdistance;
		}	
	}

    dout[myid] = minfdistance;
	
}

// Brute force - parallel 5 thread per cell method using shared memory
void ParallelBF5(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
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
    thread_per_cell_sh<<<blocks, threads, 2 * threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 5 thread per cell shared mem (GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, csize*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));


}


// Kernel for ParallelBF6: calculates the wall distance for each cell in a thread (using coalesced shared memory)
__global__ void thread_per_cell_sh_c(double * dout, int fsize, double * xc, double * yc, double * xf, double * yf)
{
	// Thread and block IDs
	int myid = threadIdx.x + blockDim.x * blockIdx.x;
    int tid  = threadIdx.x;

    // sxyf, fdistance are allocated in the kernel call: 3rd arg to <<<b, t, shmem>>>
	extern __shared__ double sdata[];
	double *sxyf = &sdata[0]; // This will contain the x and y elements interspersed, i.e. x0, y0, x1, y1, ....
	
    // Load shared mem from global mem; sxyf will contain the x and y elements interspersed, i.e. x0, y0, x1, y1, ....
	sxyf[tid*2] = xf[tid];
	sxyf[tid*2 + 1] = yf[tid];	
    __syncthreads();

	double minfdistance = 1e9; // Initialize it to a large value;
	double fdistance;

	// Loop through each face, and calculate its distance from the cell center
	for (int j=0; j<fsize; j++)
	{
//			minfdistance = min(minfdistance,sqrt(pow(xc[myid]-xf[j],2) + pow(yc[myid]-yf[j],2)));
		fdistance = sqrt(pow(xc[myid]-sxyf[j*2],2) + pow(yc[myid]-sxyf[j*2 + 1],2));
		if ( fdistance < minfdistance)
		{
			minfdistance = fdistance;
		}	
	}

    dout[myid] = minfdistance;
	
}

// Brute force - parallel 6 thread per cell method using shared memory and coalesced reads within a block
void ParallelBF6(int csize, int fsize, double * xc, double * yc, double * xf, double * yf, double * wallDist)
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
    thread_per_cell_sh_c<<<blocks, threads, 2 * threads * sizeof(double)>>>(d_wallDist, fsize, d_xc, d_yc, d_xf, d_yf);	
	timer.Stop(); // Stop timer
	printf("Brute force - parallel 6 thread per cell coalesced (GpuTimer): \t %.0f milliseconds\n",timer.Elapsed());

    // Copy result from device back to host
    checkCudaErrors(cudaMemcpy(wallDist, d_wallDist, csize*sizeof(double), cudaMemcpyDeviceToHost));

    // Clean-up memory allocated
	checkCudaErrors(cudaFree(d_xc));
    checkCudaErrors(cudaFree(d_yc));
    checkCudaErrors(cudaFree(d_xf));
    checkCudaErrors(cudaFree(d_yf));
    checkCudaErrors(cudaFree(d_wallDist));


}
