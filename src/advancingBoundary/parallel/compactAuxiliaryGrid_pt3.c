/*
 * compactAuxiliaryGrid.c
 *
 *  Created on: Nov 12, 2014
 *      Author: nathan
 */

#include <stdlib.h>
#include <stdio.h>
#include "compactAuxiliaryGrid_pt3.h"

void compactAuxiliaryGrid_pt3(struct cell_pt3 * auxCells, int numAuxCells, struct cell_pt3 * compAuxCells, double * xf, double * yf, int numFaces){

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

				index=auxCells[i].numFaces;

				// Store face and increase face count
				auxCells[i].faces[2*index] = xf[j];
				auxCells[i].faces[2*index+1] = yf[j];
				auxCells[i].numFaces++;

			}

		}
	}


	// Copy cells that contain geometry faces to compacted list
	j = 0;
	for (i=0; i<numAuxCells; i++){
		if(auxCells[i].numFaces != 0){

			// Copy contents
			compAuxCells[j].xmin = auxCells[i].xmin;
			compAuxCells[j].xmax = auxCells[i].xmax;
			compAuxCells[j].ymin = auxCells[i].ymin;
			compAuxCells[j].ymax = auxCells[i].ymax;
			compAuxCells[j].xcenter = auxCells[i].xcenter;
			compAuxCells[j].ycenter = auxCells[i].ycenter;
			compAuxCells[j].numFaces = auxCells[i].numFaces;

			for(k=0; k<2*compAuxCells[j].numFaces; k++){
				compAuxCells[j].faces[k] = auxCells[i].faces[k];
			}
			j++;
		}
	}

}
