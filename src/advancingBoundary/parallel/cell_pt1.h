/*
 * node.h
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#ifndef NODE_PT1_H_
#define NODE_PT1_H_

struct cell_pt1{

	// Cell Boundaries
	double xmin, xmax;
	double ymin, ymax;
	double xcenter, ycenter;

	// Included faces linked list
	int	numFaces;

	int storage;
	int faceIndex[20];

};



#endif /* NODE_PT1_H_ */
