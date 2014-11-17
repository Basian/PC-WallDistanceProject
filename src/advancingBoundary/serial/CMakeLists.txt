CUDA_ADD_LIBRARY( cudaTest cudaTest.cu )
#ADD_LIBRARY( ABserial OBJECT boundBox.c computeAuxiliaryGrid.c
#computeAuxiliaryGrid_t2.c compactAuxiliaryGrid.c compactAuxiliaryGrid_t2.c
#writecell.c writefaces.c ab_serial.c ab_serial_t2.c)

ADD_LIBRARY( ABserial STATIC boundBox.c computeAuxiliaryGrid.c
computeAuxiliaryGrid_t2.c compactAuxiliaryGrid.c compactAuxiliaryGrid_t2.c
writecell.c writefaces.c ab_serial.c ab_serial_t2.c)


