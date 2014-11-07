cmake_minimum_required(VERSION 3.0.1)
PROJECT(walldistancecalc Fortran CXX)

# Add preprocessor 
ADD_SUBDIRECTORY(preproc)
ADD_SUBDIRECTORY(radius)
ADD_SUBDIRECTORY(postproc)

# Add main program
ADD_EXECUTABLE( walldistcalc main.cxx )
TARGET_LINK_LIBRARIES( walldistcalc Preproc Radius Postproc)
