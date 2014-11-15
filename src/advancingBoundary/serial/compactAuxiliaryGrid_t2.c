/*
 * compactAuxiliaryGrid.c
 *
 *  Created on: Nov 12, 2014
 *      Author: nathan
 */

#include <stdlib.h>
#include <stdio.h>
#include "compactAuxiliaryGrid_t2.h"

void compactAuxiliaryGrid_t2(struct cell_t2 * auxCells, int numAuxCells, struct cell_t2 * compAuxCells, double * xf, double * yf, int numFaces){

	int i, j, k;
	int index;
	int count=0;
	int test;
	double *tempone, *temptwo;

	/* For each auxiliary cell, search through the face list and check for any faces
	 * that lie within the auxiliary cell boundary. If a face lies within the boundary,
	 * add that face as a node to the auxiliary cell linked-list of faces
	 */
	for (i=0; i<numAuxCells; i++){
		for (j=0; j<numFaces; j++){

			// If the cell contains a geometry face
			if (xf[j] < auxCells[i].xmax && xf[j] > auxCells[i].xmin && yf[j] < auxCells[i].ymax && yf[j] > auxCells[i].ymin){

				index=auxCells[i].faceNum;

				// Reallocate face storage if number of faces has exceeded current capacity
				if(index>=auxCells[i].storage){
					tempone = (double *)realloc(auxCells[i].xface, 2*auxCells[i].storage*sizeof(double));
					if(tempone == NULL){
						printf("warning, increased storage allocation failed");
					}
					else{
						auxCells[i].xface = tempone;
					}
					tempone = (double *)realloc(auxCells[i].yface, 2*auxCells[i].storage*sizeof(double));
					if(tempone == NULL){
						printf("warning, increased storage allocation failed");

					}
					else{
						auxCells[i].yface = tempone;
					}
					auxCells[i].storage *= 2;
				}

				// Store face and increase face count
				auxCells[i].xface[index] = xf[j];
				auxCells[i].yface[index] = yf[j];
				auxCells[i].faceNum++;

			}

		}
	}


	// Copy cells that contain geometry faces to compacted list
	j = 0;
	for (i=0; i<numAuxCells; i++){
		if(auxCells[i].faceNum != 0){

			// Allocate face array storage in compAuxCells
			compAuxCells[j].xface = (double *)malloc(auxCells[i].storage*sizeof(double));
			compAuxCells[j].yface = (double *)malloc(auxCells[i].storage*sizeof(double));

			// Copy contents
			compAuxCells[j].xmin = auxCells[i].xmin;
			compAuxCells[j].xmax = auxCells[i].xmax;
			compAuxCells[j].ymin = auxCells[i].ymin;
			compAuxCells[j].ymax = auxCells[i].ymax;
			compAuxCells[j].xcenter = auxCells[i].xcenter;
			compAuxCells[j].ycenter = auxCells[i].ycenter;
			compAuxCells[j].faceNum = auxCells[i].faceNum;
			compAuxCells[j].storage = auxCells[i].storage;

			for(k=0; k<compAuxCells[j].storage; k++){
				compAuxCells[j].xface[k] = auxCells[i].xface[k];
				compAuxCells[j].yface[k] = auxCells[i].yface[k];
			}
			j++;
		}
	}

}
