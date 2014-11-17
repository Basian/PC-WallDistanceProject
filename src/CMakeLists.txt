cmake_minimum_required(VERSION 2.0)
PROJECT(walldistancecalc Fortran CXX C)



find_package(CUDA)
if (CUDA_FOUND)
	include_directories(${CUDA_INCLUDE_DIRS})
endif (CUDA_FOUND)

set(CUDA_NVCC_FLAGS "-arch;compute_20" )
#set(CUDA_NVCC_FLAGS "-arch;compute_13" )
#PROJECT(walldistancecalc Fortran CXX C CUDA)


# Copy grid to build dir
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/grid/grid.x ${EXECUTABLE_OUTPUT_PATH} COPYONLY)
CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/batch.run ${EXECUTABLE_OUTPUT_PATH} COPYONLY)

# Add module directories
ADD_SUBDIRECTORY(preproc)
ADD_SUBDIRECTORY(radius)
ADD_SUBDIRECTORY(advancingBoundary)
ADD_SUBDIRECTORY(postproc)
ADD_SUBDIRECTORY(bruteforce)

# Add main program

#ADD_EXECUTABLE( walldistcalc main.cxx )
CUDA_ADD_EXECUTABLE( walldistcalc main.cxx )
TARGET_LINK_LIBRARIES( walldistcalc BruteForce Preproc Radius
ABserial ABparallel Postproc cudaTest cudaBruteForce ${CUDA_LIBRAIRES})

# Install executable
INSTALL(TARGETS walldistcalc RUNTIME DESTINATION $(CMAKE_BINARY_DIR)/bin)
