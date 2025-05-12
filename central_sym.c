#include <stdbool.h>
#include <stdckdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "vect.h"

#define MAX_POINTS 50

typedef struct {
    int data[MAX_POINTS];
    int size;
} MyArray;
bool alreadyUsedNodes[100] = {false};
int centralSymmetryTracker[100] = {-1};

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

void cs_find_best_permutation(int depth, const int numPoints, const vect *points, ring *bestVectorLength,
                              int *bestPermutation, vect currentScaledVectorAtDepth, int *VArray,
                              int *currentPermutation, MyArray *centralSymmetryList) {
    count++;
    if (depth == numPoints / 2) {
        evaluate_permutation(currentScaledVectorAtDepth, numPoints, bestVectorLength, bestPermutation,
                             currentPermutation);
        return;
    }

    for (int candidate = 0; candidate < numPoints; ++candidate) {
        int mirror = centralSymmetryList->data[candidate];
        if (alreadyUsedNodes[candidate] || alreadyUsedNodes[mirror]) continue;

        currentPermutation[depth] = candidate;
        currentPermutation[numPoints - 1 - depth] = mirror;
        
        alreadyUsedNodes[candidate] = true;
        alreadyUsedNodes[mirror] = true;

        vect nextVectScaled;
        vect_scale(&nextVectScaled, points[candidate], VArray[depth], currentScaledVectorAtDepth);
        vect_scale(&nextVectScaled, points[mirror], VArray[numPoints - 1 - depth], nextVectScaled);

        cs_find_best_permutation(depth + 1, numPoints, points, bestVectorLength, bestPermutation, nextVectScaled,
                                 VArray, currentPermutation, centralSymmetryList);

        alreadyUsedNodes[candidate] = false;
        alreadyUsedNodes[mirror] = false;
    }
}

void get_central_symmetry(const int numPoints, const vect *points, MyArray *out) {
    out->size = 0;

    if (numPoints % 2 != 0) return;

    for (int i = 0; i < numPoints; ++i) {
        for (int j = 0; j < numPoints; ++j) {
            if (i == j) continue;

            vect sum;
            vect_add(sum, points[i], points[j]);
            ring sumSquared;
            vect_norm2(sumSquared, sum);
            ring zero_ring = {0, 0};
            if (ring_comp(sumSquared, zero_ring) == 0) {
                out->data[out->size++] = j;
                break;
            }
        }
    }

    if (out->size != numPoints) {
        out->size = 0;
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

    MyArray centralSymmetryList;
    get_central_symmetry(numPoints, points, &centralSymmetryList);
    printf("Is centrally symmetric: %d\n", centralSymmetryList.size != 0);

    for (int i = 0; i < centralSymmetryList.size; ++i) {
        printf("%d ", centralSymmetryList.data[i]);
    }
    printf("\n");

    generate_zero_sum_array(numPoints, VArray);

    ring bestVectorLength = {-1, -1};
    vect currentScaledVectorAtDepth = {{0, 0}, {0, 0}, {0, 0}};

    cs_find_best_permutation(0, numPoints, points, &bestVectorLength, bestPermutation, currentScaledVectorAtDepth,
                             VArray, currentPermutation, &centralSymmetryList);

    printf("Best squared length: (%lld, %lld)\n", bestVectorLength[0], bestVectorLength[1]);
    printf("Best permutation: ");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", bestPermutation[i]);
    }
    printf("\n");

    printf("Count: %d\n", count);
}