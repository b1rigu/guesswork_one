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

int compare_norm2_desc(const void *a, const void *b) {
    const PointWithNorm *pa = a, *pb = b;
    return -ring_comp(pa->norm2, pb->norm2);  // Descending
}

int compare_coeff2_desc(const void *a, const void *b) {
    const CoeffWithSquare *ca = a, *cb = b;
    return cb->coeffSquared - ca->coeffSquared;
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
                if (newRefColumn == -1) {
                    out->data[out->size++] = candidate;
                    refColumn = newRefColumn;
                    resultDecidedNodes[candidate] = true;
                }
            }
        }

        if (newRefColumn == -1) {
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

static void vect_scale(vect *res, const vect point, int numPoints, int VArrayIdx, vect vector) {
    vect tmp;
    for (int j = 0; j < 3; j++) {
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

static long long ring_sqrt_estimate(const ring r) {
    // Convert ring = [high, low] to a 64-bit float estimate of norm²
    // Use double to safely hold all 64 bits
    double combined = (double)r[0] * (1ULL << 32) * (1ULL << 32) + (double)r[1];
    return (long long)sqrt(combined);  // Safe overestimate
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
        if ((*bestVectorLength)[0] != -1 || (*bestVectorLength)[1] != -1) {
            // vect partial_sum;
            // vect_scale(&partial_sum, points[chosen], numPoints, depth, currentScaledVectorAtDepth);

            // ring p2;
            // vect_norm2(p2, partial_sum);
            // long long P = ring_sqrt_estimate(p2);

            // long long S = 0;

            // int pointsLeftToChoose = 0;
            // int pointsIdxList[MAX_POINTS];

            // for (int g = 0; g < numPoints; ++g) {
            //     if (alreadyUsedNodes[g] || g == chosen) continue;
            //     pointsIdxList[pointsLeftToChoose] = g;
            //     pointsLeftToChoose++;
            // }

            // for (int g = 0; g < pointsLeftToChoose; ++g) {
            //     int pointIdx = pointsIdxList[g];

            //     ring v2;
            //     vect_norm2(v2, points[pointIdx]);
            //     long long V = ring_sqrt_estimate(v2);

            //     S += llabs(VArray[depth + 1 + g]) * V;
            // }

            // long long total = P + S;

            // ring bound = {total * total, 0};

            // if (ring_comp(bound, *bestVectorLength) <= 0) {
            //     continue;
            // }

            // OLD VERSION TO PRUNE
            // int pointsLeftToChoose = 0;
            // int VArrayLeftToCalcCount = 0;

            // for (int g = 0; g < numPoints; ++g) {
            //     if (g > depth) {
            //         coeffData[VArrayLeftToCalcCount].coeffIdx = g;
            //         coeffData[VArrayLeftToCalcCount].coeffSquared = VArray[g] * VArray[g];
            //         VArrayLeftToCalcCount++;
            //     }
            //     if (alreadyUsedNodes[g] || g == chosen) continue;
            //     pointData[pointsLeftToChoose].pointIndex = g;
            //     vect_norm2(pointData[pointsLeftToChoose].norm2, points[g]);
            //     pointsLeftToChoose++;
            // }

            // if (pointsLeftToChoose != VArrayLeftToCalcCount) {
            //     puts("count mismatch");
            //     exit(1);
            // }

            // qsort(&coeffData, VArrayLeftToCalcCount, sizeof(CoeffWithSquare), compare_coeff2_desc);
            // qsort(&pointData, pointsLeftToChoose, sizeof(PointWithNorm), compare_norm2_desc);

            // // printf("coeffDataSquared: ");
            // // for (int g = 0; g < VArrayLeftToCalcCount; ++g) {
            // //     printf("%d ", coeffData[g].coeffSquared);
            // // }
            // // puts("");

            // // printf("pointData: ");
            // // for (int g = 0; g < pointsLeftToChoose; ++g) {
            // //     printf("[%lld, %lld], ", pointData[g].norm2[0], pointData[g].norm2[1]);
            // // }
            // // puts("");

            // vect greedy_sum = {{0, 0}, {0, 0}, {0, 0}};

            // for (int g = 0; g < pointsLeftToChoose; ++g) {
            //     int coeffIdx = coeffData[g].coeffIdx;
            //     int pointIdx = pointData[g].pointIndex;
            //     vect_scale(&greedy_sum, points[pointIdx], numPoints, coeffIdx, greedy_sum);
            // }

            // vect total_sum;
            // vect_add(total_sum, partial_sum, greedy_sum);

            // ring total;
            // vect_norm2(total, total_sum);

            // // printf("total: [%lld, %lld]\n", total[0], total[1]);
            // // printf("current best: [%lld, %lld]\n\n", (*bestVectorLength)[0], (*bestVectorLength)[1]);

            // if (ring_comp(total, *bestVectorLength) < 0) {
            //     continue;
            // }
        }
        alreadyUsedNodes[chosen] = true;
        currentPermutation[depth] = chosen;

        // printf("best vector length: ");
        // printf("ring: [%lld, %lld]\n", (*bestVectorLength)[0], (*bestVectorLength)[1]);
        // puts("");
        // printf("Nodes to explore at depth %d :", depth);
        // for (int i = 0; i < nodes_to_explore.size; ++i) printf(" %d", nodes_to_explore.data[i]);
        // puts("");
        // printf("Chosen node: %d\n", chosen);

        vect nextVectScaled;
        vect_scale(&nextVectScaled, points[chosen], numPoints, depth, currentScaledVectorAtDepth);

        // printf("nextVectScaled: ");
        // for (int i = 0; i < 3; i++) {
        //     printf("%lld ", nextVectScaled[i][0]);
        //     printf("%lld ", nextVectScaled[i][1]);
        // }
        // puts("");

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