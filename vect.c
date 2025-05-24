#include "vect.h"
#include <stdckdint.h>
#include <stdlib.h>

void vect_add(vect res, const vect v0, const vect v1) {
    for (int i = 0; i < 3; i++) {
        ring_add(res[i], v0[i], v1[i]);
    }
}

void vect_sub(vect res, const vect v0, const vect v1) {
    for (int i = 0; i < 3; i++) {
        ring_sub(res[i], v0[i], v1[i]);
    }
}

void vect_mul(vect res, const ring x, const vect v) {
    for (int i = 0; i < 3; i++) {
        ring_mul(res[i], x, v[i]);
    }
}

void vect_norm2(ring res, const vect v)  // v = { {1, 2}, {3, 4}, {5, 6} }
{
    res[0] = 0;
    res[1] = 0;
    ring tmp;
    for (int i = 0; i < 3; i++) {
        ring_mul(tmp, v[i], v[i]);
        ring_add(res, res, tmp);
    }
}

void vect_scale(vect res, const vect point, int coeff, const vect vector) {
    vect tmp;
    for (int j = 0; j < 3; j++) {
        ring_scale(tmp[j], point[j], coeff, coeff);
    }
    vect_add(res, vector, tmp);
}