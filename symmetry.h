#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "ring.h"
#include "vect.h"

#define MAX_POINTS 50
#define MAX_SYM_GROUP_SIZE 200

void use_symmetry(int numPoints, vect *points, ring *bestVectorLength, int *bestPermutation);

#endif