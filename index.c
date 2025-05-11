#include <stdio.h>
#include <stdlib.h>

#include "symmetry.h"
#include "utils.h"

vect points[MAX_POINTS];
int bestPermutation[MAX_POINTS];
int VArray[MAX_POINTS];

int main() {
    int numPoints = nextInt();
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < 3; j++) {
            points[i][j][0] = nextLL();  // a part
            points[i][j][1] = nextLL();  // b part
        }
    }

    generate_zero_sum_array(numPoints, VArray);

    ring bestVectorLength = {-1, -1};

    use_symmetry(numPoints, points, &bestVectorLength, bestPermutation, VArray);

    printf("Best squared length: (%lld, %lld)\n", bestVectorLength[0], bestVectorLength[1]);
    printf("Best permutation: ");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", bestPermutation[i]);
    }
    printf("\n");
}