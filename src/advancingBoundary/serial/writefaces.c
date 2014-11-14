/*
 * writefaces.c
 *
 *  Created on: Nov 13, 2014
 *      Author: nathan
 */

#include "writefaces.h"
#include <stdlib.h>
#include <stdio.h>

void writefaces(struct cell * cellArray, int size){

	int i, j;
	struct face * traverse;
	FILE * facefile;
	facefile=fopen("faces_struct.txt","ab+");

	// Loop through cells
	for (i=0; i<size; i++){
		traverse = cellArray[i].root;

		while(traverse != NULL){
			fprintf(facefile, "%f %f\n",traverse->xf,traverse->yf);

			traverse = traverse->next;
		}

	}

}


