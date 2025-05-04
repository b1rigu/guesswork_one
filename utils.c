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