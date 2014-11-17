/*
 * node.h
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#ifndef NODE_T2_H_
#define NODE_T2_H_

struct cell_t3{

	// Cell Boundaries
	double xmin, xmax;
	double ymin, ymax;
	double xcenter, ycenter;

	// Included faces linked list
	int	faceNum;

	int storage;
	int * faceIndex;

};



#endif /* NODE_T2_H_ */
