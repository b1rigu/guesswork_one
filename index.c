#include <stdio.h>
#include <stdlib.h>

#include "symmetry.h"
#include "utils.h"
#include <stdckdint.h>
#include <limits.h>

#define INT128_MAX ((__int128)(((__uint128_t)1 << 127) - 1))
#define INT128_MIN (-((__int128)(((__uint128_t)1 << 127))))

vect points[MAX_POINTS];
int bestPermutation[MAX_POINTS];
int VArray[MAX_POINTS];

void print_128(__int128_t value) {
    if (value == 0) {
        putchar('0');
        return;
    }
    if (value < 0) {
        putchar('-');
        value = -value;
    }

    char buf[50];
    int i = 0;
    while (value > 0) {
        buf[i++] = '0' + (value % 10);
        value /= 10;
    }
    while (i--) putchar(buf[i]);
}

int main() {
    int numPoints = nextInt();
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < 3; j++) {
            points[i][j][0] = nextLL();  // a part
            points[i][j][1] = nextLL();  // b part
        }
    }

    printf("Number of points: %d\n", numPoints);
    printf("points:\n");
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%d ", points[i][j][0]);
            printf("%d ", points[i][j][1]);
        }
        printf(", ");
    }
    printf("\n");

    generate_zero_sum_array(numPoints, VArray);

    ring bestVectorLength = {-1, -1};

    // Currently only works with centrally symmetric shapes
    use_symmetry(numPoints, points, &bestVectorLength, bestPermutation, VArray);

    printf("Best squared length: ");
    printf("%lld ", bestVectorLength[0]);
    printf("%lld ", bestVectorLength[1]);
    printf("\n");
    printf("Best permutation: ");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", bestPermutation[i]);
    }
    printf("\n");
    puts("");

    // unsigned _BitInt(9) x = 100;
    // unsigned _BitInt(9) y = 200;
    // unsigned _BitInt(9) z = x + y;

    // printf("%u\n", (unsigned)z);  // Output: 300

    __int128_t a, b, result;

    a = INT128_MAX - 1;
    b = 1;

    result = a + b;

    print_128(result);
    printf("\n");

    

    // This will NOT compile if ckd_add doesn't support __int128_t
    if (ckd_add(&result, a, b)) {
        printf("Overflow detected\n");
    } else {
        print_128(result);
    }

    printf("Maximum bits for _BitInt: %d\n", BITINT_MAXWIDTH);
}