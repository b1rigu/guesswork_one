#include "symmetry.h"

#include <stdbool.h>
#include <stdckdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool alreadyUsedNodes[MAX_POINTS];
bool resultDecidedNodes[MAX_POINTS];
int VArray[MAX_POINTS];
int currentPermutation[MAX_POINTS];
int symmetryGroup[MAX_POINTS * MAX_SYM_GROUP_SIZE];

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

    memset(resultDecidedNodes, false, numPoints * sizeof(bool));

    int refColumn = -1;

    // Find the leftmost unused column
    for (int col = 0; col < numPoints; ++col) {
        if (!alreadyUsedNodes[col]) {
            refColumn = col;
            break;
        }
    }

    if (refColumn == -1) return;

    while (true) {
        bool foundFirstSymmetric = false;
        bool foundFirstNonSymmetric = false;
        int newRefColumn = -1;

        int refImages[MAX_SYM_GROUP_SIZE];
        for (int i = 0; i < symIndexListToSearchLength; ++i) {
            refImages[i] = SYM_IMAGE(symmetryGroup, symIndexListToSearch[i], numPoints, refColumn);
        }

        for (int candidate = 0; candidate < numPoints; ++candidate) {
            if (alreadyUsedNodes[candidate] || resultDecidedNodes[candidate]) continue;

            bool symmetryFound = false;

            for (int i = 0; i < symIndexListToSearchLength; ++i) {
                if (refImages[i] == candidate) {
                    symmetryFound = true;
                    break;
                }
            }

            if (symmetryFound) {
                if (!foundFirstSymmetric) {
                    out->data[out->size++] = candidate;
                    foundFirstSymmetric = true;
                }
                resultDecidedNodes[candidate] = true;
            } else {
                if (!foundFirstNonSymmetric) {
                    out->data[out->size++] = candidate;
                    foundFirstNonSymmetric = true;
                    newRefColumn = candidate;
                    resultDecidedNodes[candidate] = true;
                }
            }
        }

        if (foundFirstNonSymmetric) {
            refColumn = newRefColumn;
        } else {
            break;  // all remaining were symmetric → done
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
    }
}

static void vect_scale(vect *res, const vect point, int numPoints, int depth, vect currentScaledVectorAtDepth) {
    vect tmp;
    for (int j = 0; j < 3; j++) {
        if (ckd_mul(&tmp[j][0], VArray[depth], point[j][0]) || ckd_mul(&tmp[j][1], VArray[depth], point[j][1])) exit(1);
    }
    vect_add(*res, currentScaledVectorAtDepth, tmp);
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

static void find_best_permutation(int depth, const int numPoints, const int *symIndexListToSearch,
                                  const int symIndexListToSearchLength, const vect *points, ring *bestVectorLength,
                                  int *bestPermutation, vect currentScaledVectorAtDepth) {
    if (depth == numPoints) {
        evaluate_permutation(currentScaledVectorAtDepth, numPoints, bestVectorLength, bestPermutation);
        return;
    }

    MyArray nodes_to_explore;
    get_next_nodes_to_explore(numPoints, symIndexListToSearch, symIndexListToSearchLength, &nodes_to_explore);

    for (int i = 0; i < nodes_to_explore.size; ++i) {
        int chosen = nodes_to_explore.data[i];
        alreadyUsedNodes[chosen] = true;
        currentPermutation[depth] = chosen;

        vect nextVectScaled;
        vect_scale(&nextVectScaled, points[chosen], numPoints, depth, currentScaledVectorAtDepth);

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

    for (int i = 0; i < numPoints; ++i) {
        alreadyUsedNodes[i] = false;
    }

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
}