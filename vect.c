#include "vect.h"

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

void vect_dot(ring res, const vect v0, const vect v1) {
    res[0] = 0;
    res[1] = 0;
    ring temp;

    for (int i = 0; i < 3; ++i) {
        ring_mul(temp, v0[i], v1[i]);  // temp = vector1[i] * vector2[i]
        ring_add(res, res, temp);      // res += temp
    }
}

void vect_neg(vect res, const vect a)  // res = -a
{
    ring zero = {0, 0};
    for (int i = 0; i < 3; ++i) {
        ring_sub(res[i], zero, a[i]);
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

int vect_comp(const vect v0, const vect v1) {
    ring norm0;
    ring norm1;
    vect_norm2(norm0, v0);
    vect_norm2(norm1, v1);
    return ring_comp(norm0, norm1);
}