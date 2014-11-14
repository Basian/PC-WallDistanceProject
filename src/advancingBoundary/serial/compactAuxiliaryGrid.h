/*
 * compactAuxiliaryGrid.h
 *
 *  Created on: Nov 12, 2014
 *      Author: nathan
 */

#ifndef COMPACTAUXILIARYGRID_H_
#define COMPACTAUXILIARYGRID_H_

#include "cell.h"

void compactAuxiliaryGrid(struct cell * auxCells, int numCells, struct cell * compAuxCells, double * xf, double * yf, int numFaces);


#endif /* COMPACTAUXILIARYGRID_H_ */
