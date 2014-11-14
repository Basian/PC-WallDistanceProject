/*
 * compactAuxiliaryGrid.c
 *
 *  Created on: Nov 12, 2014
 *      Author: nathan
 */

#include <stdlib.h>
#include <stdio.h>
#include "compactAuxiliaryGrid.h"

void compactAuxiliaryGrid(struct cell * auxCells, int numAuxCells, struct cell * compAuxCells, double * xf, double * yf, int numFaces){

	int i;
	int j;
	int count=0;
	struct face * traverse;
	int test;

	/* For each auxiliary cell, search through the face list and check for any faces
	 * that lie within the auxiliary cell boundary. If a face lies within the boundary,
	 * add that face as a node to the auxiliary cell linked-list of faces
	 */
	for (i=0; i<numAuxCells; i++){
		traverse = auxCells[i].root;

		for (j=0; j<numFaces; j++){

			// If the cell contains a geometry face
			if (xf[j] < auxCells[i].xmax && xf[j] > auxCells[i].xmin && yf[j] < auxCells[i].ymax && yf[j] > auxCells[i].ymin){

				// Traverse faces list to find the end
				if (!traverse){

					// Allocate new face node
					auxCells[i].root = (struct face *) malloc(sizeof(struct face));
					traverse = auxCells[i].root;

					// Initialize new face node
					traverse->next = NULL;
					traverse->xf = xf[j];
					traverse->yf = yf[j];

				}
				else {
					while (traverse->next != NULL){
						traverse = traverse->next;
					}
					// Allocate new face node
					traverse->next = (struct face *) malloc(sizeof(struct face));
					traverse = traverse->next;

					// Initialize new face node
					traverse->next = NULL;
					traverse->xf = xf[j];
					traverse->yf = yf[j];
				}


			}


		}
	}


	// Copy cells that contain geometry faces to compacted list
	j = 0;
	for (i=0; i<numAuxCells; i++){
		if(auxCells[i].root){
			compAuxCells[j].xmin = auxCells[i].xmin;
			compAuxCells[j].xmax = auxCells[i].xmax;
			compAuxCells[j].ymin = auxCells[i].ymin;
			compAuxCells[j].ymax = auxCells[i].ymax;
			compAuxCells[j].xcenter = auxCells[i].xcenter;
			compAuxCells[j].ycenter = auxCells[i].ycenter;
			compAuxCells[j].root = auxCells[i].root;
			j++;
		}
	}

}
