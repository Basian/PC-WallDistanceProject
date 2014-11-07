/*
 * postproc.cxx
 *
 *  Created on: Nov 6, 2014
 *      Author: nathan
 */

#include "postproc.h"

// fortran declaration for writing plot3d field
extern"C"{
	void writevar_(int * ni, int * nj, double * var_c);
}




void postproc(int ni, int nj, double * wallDistance){

	// Write wall distance to Plot3D
	writevar_(&ni,&nj,wallDistance);


	return;
}
