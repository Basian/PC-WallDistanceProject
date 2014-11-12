/*
 * boundBox.c
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#include "boundBox.h"
#include <stdio.h>

// Function body
void boundBox(double * xf, double * yf, int size, double * xmin, double * xmax, double * ymin, double * ymax){
	int i;

	// Compute min/max
	*xmin = xf[0];
	*xmax = xf[0];
	*ymin = yf[0];
	*ymax = yf[0];

	for(i=1; i<size; i++){

		if (xf[i] < *xmin){
			*xmin = xf[i];
		}

		if (xf[i] > *xmax){
			*xmax = xf[i];
		}

		if (yf[i] < *ymin){
			*ymin = yf[i];
		}

		if (yf[i] > *ymax){
			*ymax = yf[i];
		}


	}

}

