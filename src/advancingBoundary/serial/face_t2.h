/*
 * face.h
 *
 *  Created on: Nov 13, 2014
 *      Author: nathan
 */

#ifndef FACE_H_
#define FACE_H_


// Linked-list structure for geometry faces that are included in an auxiliary cell
struct face{

	// Face center
	double xf, yf;

	// Next node
	struct face * next;

};



#endif /* FACE_H_ */
