/*
 * node.h
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#ifndef NODE_PT2_H_
#define NODE_PT2_H_

struct cell_pt2{

	// Cell Boundaries
	double xmin, xmax;
	double ymin, ymax;
	double xcenter, ycenter;

	// Included faces linked list
	int	numFaces;

	double face_x[100];
	double face_y[100];

};



#endif /* NODE_PT2_H_ */
