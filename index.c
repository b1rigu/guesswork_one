#include <stdio.h>
#include <stdlib.h>

#include "symmetry.h"
#include "utils.h"

vect points[MAX_POINTS];
int bestPermutation[MAX_POINTS];
int VArray[MAX_POINTS];

void print_scal(scal value) {
    if (value == 0) {
        putchar('0');
        return;
    }
    if (value < 0) {
        putchar('-');
        value = -value;
    }

    char buf[50];  // Enough for 128-bit integer
    int i = 0;
    while (value > 0) {
        buf[i++] = '0' + (value % 10);
        value /= 10;
    }
    while (i--) putchar(buf[i]);
}

void print_ring(const ring r) {
    printf("(");
    print_scal(r[0]);
    printf(", ");
    print_scal(r[1]);
    printf(")\n");
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

    use_symmetry(numPoints, points, &bestVectorLength, bestPermutation, VArray);

    printf("Best squared length: ");
    print_ring(bestVectorLength);
    printf("Best permutation: ");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", bestPermutation[i]);
    }
    printf("\n");
    puts("");
}