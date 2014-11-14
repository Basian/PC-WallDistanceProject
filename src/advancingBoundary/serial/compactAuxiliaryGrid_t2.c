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

	/* For each auxiliary cell, search through the face list and check for any faces
	 * that lie within the auxiliary cell boundary. If a face lies within the boundary,
	 * add that face as a node to the auxiliary cell linked-list of faces
	 */
	for (i=0; i<numAuxCells; i++){
		for (j=0; j<numFaces; j++){

			// If the cell contains a geometry face
			if (xf[j] < auxCells[i].xmax && xf[j] > auxCells[i].xmin && yf[j] < auxCells[i].ymax && yf[j] > auxCells[i].ymin){

				index=auxCells[i].faceNum;
				if(index>19){
					printf("warning! allocate more cell memory");
				}

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
			compAuxCells[j].xmin = auxCells[i].xmin;
			compAuxCells[j].xmax = auxCells[i].xmax;
			compAuxCells[j].ymin = auxCells[i].ymin;
			compAuxCells[j].ymax = auxCells[i].ymax;
			compAuxCells[j].xcenter = auxCells[i].xcenter;
			compAuxCells[j].ycenter = auxCells[i].ycenter;
			compAuxCells[j].faceNum = auxCells[i].faceNum;

			for(k=0; k<20; k++){
				compAuxCells[j].xface[k] = auxCells[i].xface[k];
				compAuxCells[j].yface[k] = auxCells[i].yface[k];
			}
			j++;
		}
	}

}
