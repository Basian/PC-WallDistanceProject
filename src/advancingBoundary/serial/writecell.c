/*
 * writecell.c
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#include "writecell.h"
#include <stdio.h>

void writecell(struct cell * cellArray, int size, int filename){
//void writecell(){

	int i;
	FILE * cellfile;

	if (filename == 0){
		cellfile=fopen("cells_init.txt","ab+");
	}
	else if (filename == 1){
		cellfile=fopen("cells_comp.txt","ab+");
	}


	for (i=0; i<size; i++){
		fprintf(cellfile, "# cell %i\n",i);
		fprintf(cellfile, "%f %f\n",cellArray[i].xmin,cellArray[i].ymin);
		fprintf(cellfile, "%f %f\n",cellArray[i].xmax,cellArray[i].ymin);
		fprintf(cellfile, "%f %f\n",cellArray[i].xmax,cellArray[i].ymax);
		fprintf(cellfile, "%f %f\n",cellArray[i].xmin,cellArray[i].ymax);
		fprintf(cellfile, "%f %f\n\n\n",cellArray[i].xmin,cellArray[i].ymin);

	}

}
