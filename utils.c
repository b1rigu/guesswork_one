#include "utils.h"

#include <stdlib.h>

int nextInt() {
    int x;
    if (scanf("%d", &x) != 1) {
        fprintf(stderr, "Error reading integer\n");
        exit(1);
    }
    return x;
}

long long nextLL() {
    long long x;
    if (scanf("%lld", &x) != 1) {
        fprintf(stderr, "Error reading long long integer\n");
        exit(1);
    }
    return x;
}

void generate_zero_sum_array(int numPoints, int *VArray) {
    int start = (numPoints - 1) * 1;  // largest positive value
    for (int i = 0; i < numPoints; ++i) {
        VArray[i] = start - 2 * i;
    }
}