
#include "cudaTest.h"
#include "cuda_utils.h"
#include <stdio.h>


__global__
void deviceFunction(int * d_array){

	int tid = threadIdx.x;

	d_array[tid] = tid;
}




void cudaTest(){

	int i;
	int * h_array;
	h_array = (int *)malloc(3*sizeof(int));
	h_array[0]=1;
	h_array[1]=1;
	h_array[2]=1;


	int * d_array;
	checkCudaErrors(cudaMalloc(&d_array,3*sizeof(int)));



	printf("Launching Three CUDA threads! :)");
	deviceFunction<<<1,3>>>(d_array);
	checkCudaErrors(cudaMemcpy(h_array, d_array, 3*sizeof(int), cudaMemcpyDeviceToHost));


	for (i=0; i<3; i++){
		printf("CUDA Thread ID: %i\n",h_array[i]);
	}

}



