#include "ring128.h"

#include <stdckdint.h>
#include <stdlib.h>
#include <stdio.h>

void ring_to_bigring(bigring out, const ring in) {
    out[0] = (bigscal)in[0];
    out[1] = (bigscal)in[1];
}

void bigring_add(bigring z, const bigring x, const bigring y) {
    if (ckd_add(&z[0], x[0], y[0]) || ckd_add(&z[1], x[1], y[1])) {
        printf("Overflow detected adding bigrings\n");
        exit(1);
    };
}

void bigring_sub(bigring z, const bigring x, const bigring y) {
    if (ckd_sub(&z[0], x[0], y[0]) || ckd_sub(&z[1], x[1], y[1])) {
        printf("Overflow detected subtracting bigrings\n");
        exit(1);
    }
}

void bigring_mul(bigring z, const bigring x, const bigring y) {
    scal tmp;
    if (ckd_mul(&z[0], x[1], y[1]) || ckd_mul(&z[0], z[0], our_k) || ckd_mul(&tmp, x[0], y[0]) ||
        ckd_add(&z[0], z[0], tmp)) {
        printf("Overflow detected multiplying bigrings\n");
        exit(1);
    };
    if (ckd_mul(&z[1], x[0], y[1]) || ckd_mul(&tmp, x[1], y[0]) || ckd_add(&z[1], z[1], tmp)) {
        printf("Overflow detected multiplying bigrings\n");
        exit(1);
    };
}

void bigring_scale(bigring z, const bigring x, int64_t scale_1, int64_t scale_2) {
    if (ckd_mul(&z[0], x[0], scale_1) || ckd_mul(&z[1], x[1], scale_2)) {
        printf("Overflow detected scaling bigrings\n");
        exit(1);
    }
}

int bigring_sign(const bigscal x) {
    if (x > 0)
        return 1;
    else if (x == 0)
        return 0;
    else
        return -1;
}

int bigring_comp(const bigring x, const bigring y) {
    bigscal tmp;
    bigring z;
    bigring_sub(z, x, y);

    int s0 = bigring_sign(z[0]);
    int s1 = bigring_sign(z[1]);

    if (s0 * s1 != -1)
        return bigring_sign(s0 + s1);
    else {
        if (ckd_mul(&tmp, z[0], z[0]) || ckd_mul(&z[1], z[1], z[1]) || ckd_mul(&z[1], our_k, z[1]) ||
            ckd_sub(&tmp, tmp, z[1])) {
            printf("Overflow detected comparing bigrings\n");
            exit(1);
        };  // tmp=z[0]*z[0]-k*z[1]*z[1]
        return s0 * bigring_sign(tmp);
    }
}