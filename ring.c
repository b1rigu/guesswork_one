#include "ring.h"

#include <stdckdint.h>
#include <stdlib.h>

void ring_add(ring z, const ring x, const ring y) {  // z[0]=x[0]+y[0] z[1]=x[1]+y[1]
    if (ckd_add(&z[0], x[0], y[0]) || ckd_add(&z[1], x[1], y[1])) exit(1);
}

void ring_sub(ring z, const ring x, const ring y) {  // z[0]=x[0]-y[0] z[1]=x[1]-y[1]
    if (ckd_sub(&z[0], x[0], y[0]) || ckd_sub(&z[1], x[1], y[1])) exit(1);
}

void ring_mul(ring z, const ring x, const ring y) {  // z[0]=x[0]*y[0]+k*x[1]*y[1] z[1]=x[0]*y[1]+y[0]*x[1]
    scal tmp;
    if (ckd_mul(&z[0], x[1], y[1]) || ckd_mul(&z[0], z[0], our_k) || ckd_mul(&tmp, x[0], y[0]) ||
        ckd_add(&z[0], z[0], tmp))
        exit(1);
    if (ckd_mul(&z[1], x[0], y[1]) || ckd_mul(&tmp, x[1], y[0]) || ckd_add(&z[1], z[1], tmp)) exit(1);
}

int ring_sign(const scal x) {
    if (x > 0)
        return 1;
    else if (x == 0)
        return 0;
    else
        return -1;
}

int ring_comp(const ring x, const ring y) {
    scal tmp;
    ring z;
    ring_sub(z, x, y);

    if (ring_sign(z[0]) * ring_sign(z[1]) != -1)
        return ring_sign(ring_sign(z[0]) + ring_sign(z[1]));
    else {
        if (ckd_mul(&tmp, z[0], z[0]) || ckd_mul(&z[1], z[1], z[1]) || ckd_mul(&z[1], our_k, z[1]) ||
            ckd_sub(&tmp, tmp, z[1]))
            exit(1);  // tmp=z[0]*z[0]-k*z[1]*z[1]
        return ring_sign(z[0]) * ring_sign(tmp);
    }
}
