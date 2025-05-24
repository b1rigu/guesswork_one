#include "ring128.h"

void ring_to_bigring(bigring out, const ring in) {
    out[0] = (bigscal)in[0];
    out[1] = (bigscal)in[1];
}

void bigring_add(bigring z, const bigring x, const bigring y) {
    z[0] = x[0] + y[0];
    z[1] = x[1] + y[1];
}
void bigring_sub(bigring z, const bigring x, const bigring y) {
    z[0] = x[0] - y[0];
    z[1] = x[1] - y[1];
}
void bigring_mul(bigring z, const bigring x, const bigring y) {
    z[0] = x[0] * y[0] + our_k * x[1] * y[1];  // your k value here
    z[1] = x[0] * y[1] + x[1] * y[0];
}

void bigring_scale(bigring z, const bigring x, int64_t scale_1, int64_t scale_2) {
    z[0] = x[0] * scale_1;
    z[1] = x[1] * scale_2;
}

int bigring_sign(const bigscal x) {
    if (x > 0)
        return 1;
    else if (x == 0)
        return 0;
    else
        return -1;
}

int bigring_comp(const bigring x,const bigring y) {
    bigring z;
    bigring_sub(z, x, y);

    int s0 = bigring_sign(z[0]);
    int s1 = bigring_sign(z[1]);

    if (s0 * s1 != -1)
        return bigring_sign(s0 + s1);
    else {
        bigscal tmp = z[0] * z[0] - z[1] * z[1] * our_k;
        return s0 * bigring_sign(tmp);
    }
}