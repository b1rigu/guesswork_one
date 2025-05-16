#include "symmetry.h"

#include <math.h>
#include <stdbool.h>
#include <stdckdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct SymmetryMap {
    int indexListToSearch[MAX_POINTS];
    int indexListSize;
    struct SymmetryMap *children[MAX_POINTS];
} SymmetryMap;

SymmetryMap allNodes[MAX_POINTS * MAX_POINTS * 1000];
int nextFreeNode = 0;

int program_subloops = 0;
int program_mainloops = 0;

int centralSymmetryTracker[100] = {-1};

ring max_norm2;
ring pointNorms[MAX_POINTS];
int coeff_sum_on_depths[MAX_POINTS] = {0};

// row 0 :  table[0*N + 0]  table[0*N + 1] … table[0*N + N-1]
// row 1 :  table[1*N + 0]  table[1*N + 1] … table[1*N + N-1]
// row 2 :  ...
#define SYM_IMAGE(table, row, N, column) ((table)[(size_t)(row) * (N) + (column)])

void precompute_point_norms(const vect *points, int numPoints, ring *pointNorms) {
    for (int i = 0; i < numPoints; ++i) {
        vect_norm2(pointNorms[i], points[i]);  // assumes vect_norm2(dest, vect)
    }
}

static void get_next_nodes_to_explore(const int numPoints, const int *symIndexListToSearch,
                                      const int symIndexListToSearchLength, MyArray *out, bool *alreadyUsedNodes,
                                      int *symmetryGroup) {
    out->size = 0;

    bool inSomeOrbit[MAX_POINTS] = {false};
    bool isFirstOrbit = true;
    for (int col = 0; col < numPoints; ++col) {
        if (!alreadyUsedNodes[col] && !inSomeOrbit[col]) {
            out->data[out->size++] = col;
            inSomeOrbit[col] = true;

            for (int i = 0; i < symIndexListToSearchLength; ++i) {
                int symIdx = symIndexListToSearch[i];
                int image = SYM_IMAGE(symmetryGroup, symIdx, numPoints, col);
                if (isFirstOrbit || image == col) {
                    inSomeOrbit[image] = true;
                }
            }
            isFirstOrbit = false;
        }
    }
};

static void get_next_symmetry_group(const int *symIndexListToSearch, int symIndexListToSearchLength, int numPoints,
                                    int chosen, MyArray *out, int *symmetryGroup) {
    out->size = 0;

    for (int i = 0; i < symIndexListToSearchLength; ++i) {
        int idx = symIndexListToSearch[i];
        int img = SYM_IMAGE(symmetryGroup, idx, numPoints, chosen);
        if (img == chosen) {
            out->data[out->size++] = idx;
        }
    }
}

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

static bool new_upper_bound(int depth, int numPoints, const vect *points, const vect currentScaledVectorAtDepth,
                            const ring *bestVectorLength, const int *VArray, const bool *alreadyUsedNodes,
                            const int chosen) {
    program_subloops++;
    // ! A = currentbestVectorLength - partial_sum_norm2 - (max_remaining_norm2 * remaining_coeff_2)
    // ! B = 4 * (max_remaining_norm2 * remaining_coeff_2) * partial_sum_norm2
    // ! C = A^2
    vect partial_sum;
    vect_scale(&partial_sum, points[chosen], VArray[depth], currentScaledVectorAtDepth);
    ring partial_sum_norm2;  // * a^2
    vect_norm2(partial_sum_norm2, partial_sum);

    // printf("remaining_coeff: %d\n", remaining_coeff);

    int remaining_coeff = coeff_sum_on_depths[depth + 1];
    scal remaining_coeff_2;
    if (ckd_mul(&remaining_coeff_2, remaining_coeff, remaining_coeff)) exit(1);

    // printf("remaining_coeff2: %d\n", remaining_coeff_2);

    ring lmax_norm2_mulby_coeff;
    if (ckd_mul(&lmax_norm2_mulby_coeff[0], max_norm2[0], remaining_coeff_2)) exit(1);
    if (ckd_mul(&lmax_norm2_mulby_coeff[1], max_norm2[1], remaining_coeff_2)) exit(1);

    // printf("lmax_norm2_mulby_coeff: (%lld, %lld)\n", lmax_norm2_mulby_coeff[0], lmax_norm2_mulby_coeff[1]);

    // * A = currentbestVectorLength - partial_sum_norm2 - (max_remaining_norm2 * remaining_coeff_2)
    ring sub_of_abc;
    ring_sub(sub_of_abc, *bestVectorLength, partial_sum_norm2);
    ring_sub(sub_of_abc, sub_of_abc, lmax_norm2_mulby_coeff);

    // printf("sub_of_abc: (%lld, %lld)\n", sub_of_abc[0], sub_of_abc[1]);

    ring zero = {0, 0};
    if (ring_comp(sub_of_abc, zero) < 0) {
        return true;
    }

    // * B = 4 * (max_remaining_norm2 * remaining_coeff_2) * partial_sum_norm2
    ring B;
    ring_mul(B, lmax_norm2_mulby_coeff, partial_sum_norm2);
    if (ckd_mul(&B[0], B[0], 4)) exit(1);
    if (ckd_mul(&B[1], B[1], 4)) exit(1);

    // printf("B: (%lld, %lld)\n", B[0], B[1]);

    // * C = A^2
    ring sub_of_abc_squared;
    ring_mul(sub_of_abc_squared, sub_of_abc, sub_of_abc);

    // printf("sub_of_abc_squared: (%lld, %lld)\n", sub_of_abc_squared[0], sub_of_abc_squared[1]);

    return ring_comp(sub_of_abc_squared, B) < 0;
}

static bool calculate_upper_bound(int depth, int numPoints, const vect *points, const vect currentScaledVectorAtDepth,
                                  const ring *bestVectorLength, const int *VArray, const bool *alreadyUsedNodes,
                                  const int chosen, const MyArray *centralSymmetryList) {
    program_subloops++;
    int mirror = centralSymmetryList->data[chosen];
    vect partial_sum;
    vect_scale(&partial_sum, points[chosen], VArray[depth], currentScaledVectorAtDepth);
    vect_scale(&partial_sum, points[mirror], VArray[numPoints - 1 - depth], partial_sum);

    vect max_comp;
    for (int i = 0; i < 3; ++i) {
        program_subloops++;
        max_comp[i][0] = partial_sum[i][0];
        max_comp[i][1] = partial_sum[i][1];

        scal abs_sum0 = 0;
        scal abs_sum1 = 0;

        scal temp0, temp1;

        int next_coeff_idx = depth + 1;

        for (int g = 0; g < numPoints; ++g) {
            program_subloops++;
            if (alreadyUsedNodes[g] || chosen == g || mirror == g) continue;

            int mirror_g = centralSymmetryList->data[g];

            if (next_coeff_idx >= numPoints / 2) break;

            int coeff1 = VArray[next_coeff_idx];
            int coeff2 = VArray[numPoints - 1 - next_coeff_idx];
            next_coeff_idx++;

            scal prod0a, prod1a, prod0b, prod1b;
            if (ckd_mul(&prod0a, coeff1, points[g][i][0]) || ckd_mul(&prod1a, coeff1, points[g][i][1]) ||
                ckd_mul(&prod0b, coeff2, points[mirror_g][i][0]) || ckd_mul(&prod1b, coeff2, points[mirror_g][i][1])) {
                exit(1);
            }

            // abs_sum0 += (prod0a > 0 ? prod0a : -prod0a) + (prod0b > 0 ? prod0b : -prod0b);
            // abs_sum1 += (prod1a > 0 ? prod1a : -prod1a) + (prod1b > 0 ? prod1b : -prod1b);

            if (ckd_add(&temp0, (prod0a > 0 ? prod0a : -prod0a), (prod0b > 0 ? prod0b : -prod0b))) exit(1);
            if (ckd_add(&abs_sum0, abs_sum0, temp0)) exit(1);

            if (ckd_add(&temp1, (prod1a > 0 ? prod1a : -prod1a), (prod1b > 0 ? prod1b : -prod1b))) exit(1);
            if (ckd_add(&abs_sum1, abs_sum1, temp1)) exit(1);
        }

        // max_comp[i][0] += abs_sum0;
        // max_comp[i][1] += abs_sum1;
        if (ckd_add(&max_comp[i][0], max_comp[i][0], abs_sum0)) exit(1);
        if (ckd_add(&max_comp[i][1], max_comp[i][1], abs_sum1)) exit(1);
    }

    ring comp_norm2;
    vect_norm2(comp_norm2, max_comp);

    return ring_comp(comp_norm2, *bestVectorLength) > 0;
}

static void find_best_permutation(int depth, const int numPoints, const vect *points, ring *bestVectorLength,
                                  int *bestPermutation, vect currentScaledVectorAtDepth, int *VArray,
                                  int *currentPermutation, bool *alreadyUsedNodes, SymmetryMap *root,
                                  MyArray *centralSymmetryList) {
    program_mainloops++;
    if (depth == numPoints / 2) {
        evaluate_permutation(currentScaledVectorAtDepth, numPoints, bestVectorLength, bestPermutation,
                             currentPermutation);
        return;
    }

    for (int candidate = 0; candidate < numPoints; ++candidate) {
        int mirror = centralSymmetryList->data[candidate];
        if (alreadyUsedNodes[candidate] || alreadyUsedNodes[mirror]) continue;
        if (root != NULL && root->indexListSize > 0 && root->children[candidate] == NULL) continue;
        if ((*bestVectorLength)[0] != -1 || (*bestVectorLength)[1] != -1) {
            // if (!new_upper_bound(depth, numPoints, points, currentScaledVectorAtDepth, bestVectorLength, VArray,
            //                      alreadyUsedNodes, candidate)) {
            //     continue;  // Skip this branch if upper bound is not promising
            // }
            if (!calculate_upper_bound(depth, numPoints, points, currentScaledVectorAtDepth, bestVectorLength, VArray,
                                       alreadyUsedNodes, candidate, centralSymmetryList)) {
                continue;  // Skip this branch if upper bound is not promising
            }
        }

        currentPermutation[depth] = candidate;
        currentPermutation[numPoints - 1 - depth] = mirror;

        alreadyUsedNodes[candidate] = true;
        alreadyUsedNodes[mirror] = true;

        vect nextVectScaled;
        vect_scale(&nextVectScaled, points[candidate], VArray[depth], currentScaledVectorAtDepth);
        // vect_scale(&nextVectScaled, points[mirror], VArray[numPoints - 1 - depth], nextVectScaled);

        SymmetryMap *nextChild = (root != NULL && root->indexListSize > 0) ? root->children[candidate] : NULL;

        find_best_permutation(depth + 1, numPoints, points, bestVectorLength, bestPermutation, nextVectScaled, VArray,
                              currentPermutation, alreadyUsedNodes, nextChild, centralSymmetryList);

        alreadyUsedNodes[candidate] = false;
        alreadyUsedNodes[mirror] = false;
    }
}

SymmetryMap *getNewNode() {
    if (nextFreeNode >= MAX_POINTS * MAX_POINTS * 1000) {
        // Error: out of nodes
        return NULL;
    }

    SymmetryMap *newNode = &allNodes[nextFreeNode++];
    newNode->indexListSize = 0;
    for (int i = 0; i < MAX_POINTS; i++) {
        newNode->children[i] = NULL;
    }
    return newNode;
}

static void prepare_symmetry_graph(const int numPoints, const int *symIndexListToSearch,
                                   const int symIndexListToSearchLength, bool *alreadyUsedNodes, int *symmetryGroup,
                                   SymmetryMap *root) {
    program_subloops++;
    if (symIndexListToSearchLength == 1) {
        return;
    }

    MyArray nodes_to_explore;
    get_next_nodes_to_explore(numPoints, symIndexListToSearch, symIndexListToSearchLength, &nodes_to_explore,
                              alreadyUsedNodes, symmetryGroup);

    for (int i = 0; i < nodes_to_explore.size; ++i) {
        int nextChosen = nodes_to_explore.data[i];
        alreadyUsedNodes[nextChosen] = true;

        MyArray nextSymIndexList;
        get_next_symmetry_group(symIndexListToSearch, symIndexListToSearchLength, numPoints, nextChosen,
                                &nextSymIndexList, symmetryGroup);

        SymmetryMap *child = getNewNode();
        if (child == NULL) {
            alreadyUsedNodes[nextChosen] = false;
            return;
        };

        root->children[nextChosen] = child;
        root->indexListToSearch[root->indexListSize++] = nextChosen;

        prepare_symmetry_graph(numPoints, nextSymIndexList.data, nextSymIndexList.size, alreadyUsedNodes, symmetryGroup,
                               child);

        alreadyUsedNodes[nextChosen] = false;
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

void use_symmetry(int numPoints, vect *points, ring *bestVectorLength, int *bestPermutation, int *VArray) {
    bool alreadyUsedNodes[MAX_POINTS] = {false};
    int currentPermutation[MAX_POINTS] = {0};
    int symmetryGroup[MAX_POINTS * MAX_SYM_GROUP_SIZE];
    precompute_point_norms(points, numPoints, pointNorms);

    // printf("point norms precomputed\n");
    // for (int i = 0; i < numPoints; ++i) {
    //     printf("%d ", pointNorms[i][0]);
    //     printf("%d ", pointNorms[i][1]);
    //     printf("\n");
    // }

    int total = 0;
    for (int i = numPoints / 2 - 1; i >= 0; --i) {
        total += VArray[i];
        coeff_sum_on_depths[i] = total;
    }

    printf("varray: ");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", VArray[i]);
    }
    printf("\n");

    printf("coeff sum on depths\n");
    for (int i = 0; i < numPoints / 2; ++i) {
        printf("%d ", coeff_sum_on_depths[i]);
    }
    printf("\n");

    bool first_iter = true;
    for (int i = 0; i < numPoints; ++i) {
        program_subloops++;
        if (first_iter) {
            max_norm2[0] = pointNorms[i][0];
            max_norm2[1] = pointNorms[i][1];
            first_iter = false;
            continue;
        }
        if (ring_comp(pointNorms[i], max_norm2) > 0) {
            max_norm2[0] = pointNorms[i][0];
            max_norm2[1] = pointNorms[i][1];
        }
    }

    int symmetryGroupSize = nextInt();
    for (int g = 0; g < symmetryGroupSize; g++) {
        for (int i = 0; i < numPoints; i++) {
            symmetryGroup[g * numPoints + i] = nextInt();
        }
    }
    int symIndexListToSearch[MAX_SYM_GROUP_SIZE];
    for (int g = 0; g < symmetryGroupSize; ++g) {
        symIndexListToSearch[g] = g;
    }

    MyArray centralSymmetryList;
    get_central_symmetry(numPoints, points, &centralSymmetryList);
    printf("Is centrally symmetric: %d\n", centralSymmetryList.size != 0);

    vect currentScaledVectorAtDepth = {{0, 0}, {0, 0}, {0, 0}};

    SymmetryMap *root = getNewNode();
    prepare_symmetry_graph(numPoints, symIndexListToSearch, symmetryGroupSize, alreadyUsedNodes, symmetryGroup, root);

    printf("symmetry map prep done\n");
    printf("Symmetry prep loops: %d\n", program_subloops);

    memset(alreadyUsedNodes, false, sizeof(alreadyUsedNodes));

    clock_t begin = clock();
    find_best_permutation(0, numPoints, points, bestVectorLength, bestPermutation, currentScaledVectorAtDepth, VArray,
                          currentPermutation, alreadyUsedNodes, root, &centralSymmetryList);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    (*bestVectorLength)[0] = (*bestVectorLength)[0] * 4;
    (*bestVectorLength)[1] = (*bestVectorLength)[1] * 4;

    printf("Total Runtime Count subloops: %d\n", program_subloops);
    printf("Total Runtime Count mainloops: %d\n", program_mainloops);
    printf("Time spent: %f\n", time_spent);
}