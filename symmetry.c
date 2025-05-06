#include "symmetry.h"

#include <math.h>
#include <stdbool.h>
#include <stdckdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int pointIndex;
    ring norm2;
} PointWithNorm;

typedef struct {
    int coeffIdx;      // original signed
    int coeffSquared;  // coeff * coeff
} CoeffWithSquare;

PointWithNorm pointData[MAX_POINTS];
CoeffWithSquare coeffData[MAX_POINTS];
bool alreadyUsedNodes[MAX_POINTS] = {false};
bool resultDecidedNodes[MAX_POINTS] = {false};
int VArray[MAX_POINTS] = {0};
int currentPermutation[MAX_POINTS] = {0};
int symmetryGroup[MAX_POINTS * MAX_SYM_GROUP_SIZE] = {0};
int program_subloops = 0;
int program_mainloops = 0;

// row 0 :  table[0*N + 0]  table[0*N + 1] … table[0*N + N-1]
// row 1 :  table[1*N + 0]  table[1*N + 1] … table[1*N + N-1]
// row 2 :  ...
#define SYM_IMAGE(table, row, N, column) ((table)[(size_t)(row) * (N) + (column)])

static void generate_zero_sum_array(int numPoints, int *VArray) {
    int start = (numPoints - 1) * 1;  // largest positive value
    for (int i = 0; i < numPoints; ++i) {
        VArray[i] = start - 2 * i;
    }
}

static void get_next_nodes_to_explore(const int numPoints, const int *symIndexListToSearch,
                                      const int symIndexListToSearchLength, MyArray *out) {
    out->size = 0;

    bool inSomeOrbit[MAX_POINTS] = {false};
    bool isFirstOrbit = true;
    for (int col = 0; col < numPoints; ++col) {
        if (!alreadyUsedNodes[col] && !inSomeOrbit[col]) {
            out->data[out->size++] = col;
            inSomeOrbit[col] = true;

            for (int i = 0; i < symIndexListToSearchLength; ++i) {
                program_subloops++;
                int symIdx = symIndexListToSearch[i];
                int image = SYM_IMAGE(symmetryGroup, symIdx, numPoints, col);
                if (isFirstOrbit || image == col) {
                    inSomeOrbit[image] = true;
                }
            }
            isFirstOrbit = false;
            program_subloops++;
        }
    }
};

static void get_next_symmetry_group(const int *symIndexListToSearch, int symIndexListToSearchLength, int numPoints,
                                    int chosen, MyArray *out) {
    out->size = 0;

    for (int i = 0; i < symIndexListToSearchLength; ++i) {
        int idx = symIndexListToSearch[i];
        int img = SYM_IMAGE(symmetryGroup, idx, numPoints, chosen);
        if (img == chosen) {
            out->data[out->size++] = idx;
        }
        program_subloops++;
    }
}

static void vect_scale(vect *res, const vect point, int VArrayIdx, const vect vector) {
    vect tmp;
    for (int j = 0; j < 3; j++) {
        program_subloops++;
        if (ckd_mul(&tmp[j][0], VArray[VArrayIdx], point[j][0]) || ckd_mul(&tmp[j][1], VArray[VArrayIdx], point[j][1]))
            exit(1);
    }
    vect_add(*res, vector, tmp);
}

static void evaluate_permutation(vect currentScaledVector, int numPoints, ring *bestVectorLength,
                                 int *bestPermutation) {
    ring currentSquared;
    vect_norm2(currentSquared, currentScaledVector);
    if (ring_comp(currentSquared, *bestVectorLength) > 0) {
        (*bestVectorLength)[0] = currentSquared[0];
        (*bestVectorLength)[1] = currentSquared[1];
        memcpy(bestPermutation, currentPermutation, numPoints * sizeof(*currentPermutation));
    }
}

int compare_norm2_desc(const void *a, const void *b) {
    const PointWithNorm *pa = a, *pb = b;
    return -ring_comp(pa->norm2, pb->norm2);  // Descending
}

int compare_coeff2_desc(const void *a, const void *b) {
    const CoeffWithSquare *ca = a, *cb = b;
    return cb->coeffSquared - ca->coeffSquared;
}

static bool calculate_upper_bound(int depth, const int numPoints, const vect *points,
                                  const vect currentScaledVectorAtDepth, const ring *bestVectorLength,
                                  const int chosen) {
    vect partial_sum;
    vect_scale(&partial_sum, points[chosen], depth, currentScaledVectorAtDepth);

    int pointsLeftToChoose = 0;
    int VArrayLeftToCalcCount = 0;

    for (int g = 0; g < numPoints; ++g) {
        program_subloops++;
        if (g > depth) {
            coeffData[VArrayLeftToCalcCount].coeffIdx = g;
            coeffData[VArrayLeftToCalcCount].coeffSquared = VArray[g] * VArray[g];
            VArrayLeftToCalcCount++;
        }
        if (alreadyUsedNodes[g] || g == chosen) continue;
        pointData[pointsLeftToChoose].pointIndex = g;
        vect_norm2(pointData[pointsLeftToChoose].norm2, points[g]);
        pointsLeftToChoose++;
    }

    qsort(&coeffData, VArrayLeftToCalcCount, sizeof(CoeffWithSquare), compare_coeff2_desc);
    qsort(&pointData, pointsLeftToChoose, sizeof(PointWithNorm), compare_norm2_desc);

    vect max_comp = {{0, 0}, {0, 0}, {0, 0}};
    for (int i = 0; i < 3; i++) {
        program_subloops++;
        max_comp[i][0] = partial_sum[i][0];
        max_comp[i][1] = partial_sum[i][1];

        // Calculate sum of absolute contributions for each component
        scal abs_contrib_0 = 0;
        scal abs_contrib_1 = 0;

        for (int g = 0; g < pointsLeftToChoose; ++g) {
            program_subloops++;
            int coeffPos = coeffData[g].coeffIdx;
            int pointIdx = pointData[g].pointIndex;

            scal prod_0, prod_1;
            if (ckd_mul(&prod_0, VArray[coeffPos], points[pointIdx][i][0]) ||
                ckd_mul(&prod_1, VArray[coeffPos], points[pointIdx][i][1]))
                exit(1);

            abs_contrib_0 += (prod_0 > 0) ? prod_0 : -prod_0;
            abs_contrib_1 += (prod_1 > 0) ? prod_1 : -prod_1;
        }

        // Add to the component-wise upper bound
        max_comp[i][0] += abs_contrib_0;
        max_comp[i][1] += abs_contrib_1;
    }

    ring comp_norm;
    vect_norm2(comp_norm, max_comp);

    comp_norm[0] += 1;

    return ring_comp(comp_norm, *bestVectorLength) <= 0;
}

static void find_best_permutation(int depth, const int numPoints, const int *symIndexListToSearch,
                                  const int symIndexListToSearchLength, const vect *points, ring *bestVectorLength,
                                  int *bestPermutation, vect currentScaledVectorAtDepth) {
    program_mainloops++;
    if (depth == numPoints) {
        evaluate_permutation(currentScaledVectorAtDepth, numPoints, bestVectorLength, bestPermutation);
        return;
    }

    MyArray nodes_to_explore;
    get_next_nodes_to_explore(numPoints, symIndexListToSearch, symIndexListToSearchLength, &nodes_to_explore);

    for (int i = 0; i < nodes_to_explore.size; ++i) {
        int chosen = nodes_to_explore.data[i];
        if ((*bestVectorLength)[0] != -1 || (*bestVectorLength)[1] != -1) {
            if (calculate_upper_bound(depth, numPoints, points, currentScaledVectorAtDepth, bestVectorLength,
            chosen)) {
                continue;  // Skip this branch if upper bound is not promising
            }
        }
        alreadyUsedNodes[chosen] = true;
        currentPermutation[depth] = chosen;

        vect nextVectScaled;
        vect_scale(&nextVectScaled, points[chosen], depth, currentScaledVectorAtDepth);

        MyArray nextSymIndexList;
        get_next_symmetry_group(symIndexListToSearch, symIndexListToSearchLength, numPoints, chosen, &nextSymIndexList);

        find_best_permutation(depth + 1, numPoints, nextSymIndexList.data, nextSymIndexList.size, points,
                              bestVectorLength, bestPermutation, nextVectScaled);

        alreadyUsedNodes[chosen] = false;
    }
}

void use_symmetry(int numPoints, vect *points, ring *bestVectorLength, int *bestPermutation) {
    int symmetryGroupSize = nextInt();
    for (int g = 0; g < symmetryGroupSize; g++) {
        for (int i = 0; i < numPoints; i++) {
            symmetryGroup[g * numPoints + i] = nextInt();
        }
    }
    // Print symmetry group permutations
    printf("Symmetry group list:\n");
    for (int g = 0; g < symmetryGroupSize; g++) {
        for (int i = 0; i < numPoints; i++) {
            printf("%d ", symmetryGroup[g * numPoints + i]);
        }
        printf(", ");
    }
    printf("\n");

    // Print N points
    printf("N points:\n");
    for (int i = 0; i < numPoints; i++) {
        printf("Point %d: ", i);
        for (int j = 0; j < 3; j++) {
            printf("%lld ", points[i][j][0]);
            printf("%lld ", points[i][j][1]);
        }
        printf("\n");
    }
    printf("\n");

    int symIndexListToSearch[MAX_SYM_GROUP_SIZE];
    for (int g = 0; g < symmetryGroupSize; ++g) {
        symIndexListToSearch[g] = g;
    }
    generate_zero_sum_array(numPoints, VArray);
    printf("VArray: [");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", VArray[i]);
    }
    puts("]");

    vect currentScaledVectorAtDepth = {{0, 0}, {0, 0}, {0, 0}};
    find_best_permutation(0, numPoints, symIndexListToSearch, symmetryGroupSize, points, bestVectorLength,
                          bestPermutation, currentScaledVectorAtDepth);

    printf("Total Runtime Count subloops: %d\n", program_subloops);
    printf("Total Runtime Count mainloops: %d\n", program_mainloops);
}