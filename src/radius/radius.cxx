/*
 * radius.cxx
 *
 *  Created on: Nov 7, 2014
 *      Author: nathan
 */

#include <math.h>
#include "radius.h"

void radius(int size, double * xc, double * yc, double * r){

	for (int i=0; i<size; i++){
		r[i] = sqrt(pow(xc[i],2) + pow(yc[i],2));
	}

}


