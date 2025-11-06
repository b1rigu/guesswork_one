#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

typedef __uint64_t bitmask_t;

#define IS_USED(mask, i) (((mask) >> (i)) & 1)
#define SET_USED(mask, i) ((mask) |= ((bitmask_t)1 << (i)))
#define UNSET_USED(mask, i) ((mask) &= ~((bitmask_t)1 << (i)))

// Function to safely read an integer from input
int nextInt();

// Function to safely read a long long integer from input
long long nextLL();

#endif
