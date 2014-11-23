/*
 * node.h
 *
 *  Created on: Nov 11, 2014
 *      Author: nathan
 */

#ifndef NODE_PT3_H_
#define NODE_PT3_H_

struct cell_pt3{

	// Cell Boundaries
	double xmin, xmax;
	double ymin, ymax;
	double xcenter, ycenter;

	// Included faces linked list
	int	numFaces;

	double faces[200];

};



#endif /* NODE_PT3_H_ */
