#include <stdbool.h>
#include <stdckdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "vect.h"

bool alreadyUsedNodes[100] = {false};

static void evaluate_permutation(vect currentScaledVector, int numPoints, ring *bestVectorLength, int *bestPermutation,
                                 int *currentPermutation) {
    ring currentSquared;
    vect_norm2(currentSquared, currentScaledVector);
    if (ring_comp(currentSquared, *bestVectorLength) > 0) {
        (*bestVectorLength)[0] = currentSquared[0];
        (*bestVectorLength)[1] = currentSquared[1];
        memcpy(bestPermutation, currentPermutation, numPoints * sizeof(*currentPermutation));
    }
}

int count = 0;

void brute_force_find_best_permutation(int depth, const int numPoints, const vect *points, ring *bestVectorLength,
                                       int *bestPermutation, vect currentScaledVectorAtDepth, int *VArray,
                                       int *currentPermutation) {
    count++;
    if (depth == numPoints) {
        evaluate_permutation(currentScaledVectorAtDepth, numPoints, bestVectorLength, bestPermutation,
                             currentPermutation);
        return;
    }

    for (int candidate = 0; candidate < numPoints; ++candidate) {
        if (alreadyUsedNodes[candidate]) continue;
        alreadyUsedNodes[candidate] = true;
        currentPermutation[depth] = candidate;

        vect nextVectScaled;
        vect_scale(&nextVectScaled, points[candidate], VArray[depth], currentScaledVectorAtDepth);

        brute_force_find_best_permutation(depth + 1, numPoints, points, bestVectorLength, bestPermutation,
                                          nextVectScaled, VArray, currentPermutation);

        alreadyUsedNodes[candidate] = false;
    }
}

int main() {
    int currentPermutation[100] = {0};
    int VArray[100];
    int bestPermutation[100];
    vect points[100];

    int numPoints = nextInt();
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < 3; j++) {
            points[i][j][0] = nextLL();  // a part
            points[i][j][1] = nextLL();  // b part
        }
    }

    generate_zero_sum_array(numPoints, VArray);

    ring bestVectorLength = {-1, -1};
    vect currentScaledVectorAtDepth = {{0, 0}, {0, 0}, {0, 0}};

    brute_force_find_best_permutation(0, numPoints, points, &bestVectorLength, bestPermutation,
                                      currentScaledVectorAtDepth, VArray, currentPermutation);

    printf("Best squared length: (%lld, %lld)\n", bestVectorLength[0], bestVectorLength[1]);
    printf("Best permutation: ");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", bestPermutation[i]);
    }
    printf("\n");

    printf("Count: %d\n", count);
}