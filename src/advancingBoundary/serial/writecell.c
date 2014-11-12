/*
 * writecell.c
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#include "writecell.h"
#include <stdio.h>

void writecell(struct cell * cellArray, int size){
//void writecell(){

	int i;

	FILE * cellfile=fopen("cells.txt","ab+");

	for (i=0; i<size; i++){
		fprintf(cellfile, "# cell %i\n",i);
		fprintf(cellfile, "%f %f\n",cellArray[i].xmin,cellArray[i].ymin);
		fprintf(cellfile, "%f %f\n",cellArray[i].xmax,cellArray[i].ymin);
		fprintf(cellfile, "%f %f\n",cellArray[i].xmax,cellArray[i].ymax);
		fprintf(cellfile, "%f %f\n",cellArray[i].xmin,cellArray[i].ymax);
		fprintf(cellfile, "%f %f\n",cellArray[i].xmin,cellArray[i].ymin);

	}






}
