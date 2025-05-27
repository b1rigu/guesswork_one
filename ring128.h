#ifndef RING128_H
#define RING128_H

#include <stdint.h>
#include "ring.h"

typedef __int128_t bigscal;
typedef bigscal bigring[2];

void ring_to_bigring(bigring out, const ring in);

void bigring_add(bigring z, const bigring x, const bigring y);

void bigring_sub(bigring z, const bigring x, const bigring y);

void bigring_mul(bigring z, const bigring x, const bigring y);

void bigring_scale(bigring z, const bigring x, int64_t scale_1, int64_t scale_2);

int bigring_sign(const bigscal x);

int bigring_comp(const bigring x,const bigring y);

#endif
