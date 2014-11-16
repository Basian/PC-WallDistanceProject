ADD_LIBRARY( BruteForce STATIC bruteforce.cxx )
CUDA_ADD_LIBRARY( cudaBruteForce bruteforce_parallel.cu )
