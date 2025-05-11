#ifndef VECT_H
#define VECT_H

#include "ring.h"

typedef ring vect[3];

void vect_add(vect res, const vect v0, const vect v1);

void vect_sub(vect res, const vect v0, const vect v1);

void vect_mul(vect res, const ring x, const vect v);

void vect_dot(ring res, const vect v0, const vect v1);

void vect_neg(vect res, const vect a);

int vect_comp(const vect v0, const vect v1);

void vect_norm2(ring res, const vect v);

void vect_scale(vect *res, const vect point, int coeff, const vect vector);

#endif
