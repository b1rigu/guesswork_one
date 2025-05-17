#ifndef RING_H
#define RING_H

#define our_k 5

typedef __int128_t scal;

typedef scal ring[2];

void ring_add(ring z, const ring x, const ring y);

void ring_sub(ring z, const ring x, const ring y);

void ring_mul(ring z, const ring x, const ring y);

void ring_scale(ring z, const ring x, int scale_1, int scale_2);

int ring_sign(const scal x);

int ring_comp(const ring x,const ring y);

#endif
