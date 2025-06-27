#include <stdckdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "symmetry.h"
#include "utils.h"

int main() {
    vect points[MAX_POINTS];
    int bestPermutation[MAX_POINTS];

    int numPoints = nextInt();
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < 3; j++) {
            points[i][j][0] = nextLL();  // a part
            points[i][j][1] = nextLL();  // b part
        }
    }

    ring bestVectorLength = {-1, -1};

    // Currently only works with centrally symmetric shapes
    use_symmetry(numPoints, points, &bestVectorLength, bestPermutation);

    printf("Best squared length: ");
    printf("%lld ", bestVectorLength[0]);
    printf("%lld ", bestVectorLength[1]);
    printf("\n");
    printf("Best permutation: ");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", bestPermutation[i]);
    }
    printf("\n");
}