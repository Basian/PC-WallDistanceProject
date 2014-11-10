cmake_minimum_required(VERSION 2.0)
PROJECT(walldistancecalc Fortran CXX)

# Copy grid to build dir
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/grid/grid.x ${EXECUTABLE_OUTPUT_PATH} COPYONLY)

# Add module directories
ADD_SUBDIRECTORY(preproc)
ADD_SUBDIRECTORY(radius)
ADD_SUBDIRECTORY(postproc)

# Add main program
ADD_EXECUTABLE( walldistcalc main.cxx )
TARGET_LINK_LIBRARIES( walldistcalc Preproc Radius Postproc)

# Install executable
INSTALL(TARGETS walldistcalc RUNTIME DESTINATION $(CMAKE_BINARY_DIR)/bin)
