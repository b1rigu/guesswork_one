#include "symmetry.h"

#include <math.h>
#include <stdbool.h>
#include <stdckdint.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct SymmetryMap {
    signed char children_index_list[MAX_POINTS];
    unsigned char children_count;
    struct SymmetryMap **children;
} SymmetryMap;

typedef __uint64_t bitmask_t;

static double sqrt_our_k = 0.0;
int best_available_index = 0;

#define IS_USED(mask, i) (((mask) >> (i)) & 1)
#define SET_USED(mask, i) ((mask) |= ((bitmask_t)1 << (i)))
#define UNSET_USED(mask, i) ((mask) &= ~((bitmask_t)1 << (i)))

// row 0 :  table[0*N + 0]  table[0*N + 1] … table[0*N + N-1]
// row 1 :  table[1*N + 0]  table[1*N + 1] … table[1*N + N-1]
// row 2 :  ...
#define SYM_IMAGE(table, row, N, column) ((table)[(size_t)(row) * (N) + (column)])

SymmetryMap *createNode() {
    SymmetryMap *node = malloc(sizeof(SymmetryMap));
    if (!node) {
        perror("malloc failed");
        exit(1);
    }

    for (int i = 0; i < MAX_POINTS; i++) {
        node->children_index_list[i] = -1;
    }

    node->children_count = 0;
    node->children = NULL;

    return node;
}

void addChild(SymmetryMap *node, unsigned char index, SymmetryMap *child) {
    signed char pos = node->children_index_list[index];

    if (pos != -1) {
        // Overwrite existing child
        node->children[pos] = child;
        return;
    }

    // Add new child
    unsigned char newIndex = node->children_count;

    node->children = realloc(node->children, sizeof(SymmetryMap *) * (newIndex + 1));
    if (!node->children) {
        perror("realloc failed");
        exit(1);
    }

    node->children[newIndex] = child;
    node->children_index_list[index] = newIndex;
    node->children_count++;
}

static SymmetryMap *getChild(SymmetryMap *node, unsigned char index) {
    if (node == NULL) return NULL;
    signed char pos = node->children_index_list[index];
    if (pos == -1) return NULL;
    return node->children[pos];
}

void freeSymmetryMap(SymmetryMap *node) {
    if (!node) return;

    for (int i = 0; i < node->children_count; i++) {
        freeSymmetryMap(node->children[i]);
    }

    free(node->children);
    free(node);
}

static void precompute_point_norms(const vect *points, int numPoints, ring *pointNorms) {
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

static bool cheap_upper_bound_estimate(int depth, const vect nextScaledVect, const ring *bestVectorLength,
                                       int32_t *coeff_sum_on_depths, const ring *pointNorms) {
    ring partial_sum_norm2;
    vect_norm2(partial_sum_norm2, nextScaledVect);

    const ring *current_max_norm2 = &pointNorms[best_available_index];

    double remaining_coeff = (double)coeff_sum_on_depths[depth + 1];

    double best_sq = (double)(*bestVectorLength)[0] + (double)(*bestVectorLength)[1] * sqrt_our_k;
    double partial_sq = (double)partial_sum_norm2[0] + (double)partial_sum_norm2[1] * sqrt_our_k;
    double max_sq = (double)(*current_max_norm2)[0] + (double)(*current_max_norm2)[1] * sqrt_our_k;

    double total_est = sqrt(partial_sq) + remaining_coeff * sqrt(max_sq);
    return total_est > sqrt(best_sq);
}

static bool should_check_branch_ub(int depth, const vect nextScaledVect, const ring *bestVectorLength,
                                   const ring *pointNorms, int32_t *coeff_sum_on_depths) {
    // ! A = currentbestVectorLength - partial_sum_norm2 - (max_remaining_norm2 * remaining_coeff_2)
    // ! B = 4 * (max_remaining_norm2 * remaining_coeff_2) * partial_sum_norm2
    // ! C = A^2
    ring partial_sum_norm2;  // * a^2
    vect_norm2(partial_sum_norm2, nextScaledVect);
    bigring partial_sum_norm2_big;
    ring_to_bigring(partial_sum_norm2_big, partial_sum_norm2);

    // printf("remaining_coeff: %d\n", remaining_coeff);

    int64_t remaining_coeff = coeff_sum_on_depths[depth + 1];
    int64_t remaining_coeff_2 = remaining_coeff * remaining_coeff;

    // printf("remaining_coeff2: %d\n", remaining_coeff_2);

    const ring *current_max_norm2 = &pointNorms[best_available_index];

    bigring lmax_norm2_mulby_coeff;
    bigring current_max_norm2_big;
    ring_to_bigring(current_max_norm2_big, *current_max_norm2);
    bigring_scale(lmax_norm2_mulby_coeff, current_max_norm2_big, remaining_coeff_2, remaining_coeff_2);

    // printf("lmax_norm2_mulby_coeff: (%lld, %lld)\n", lmax_norm2_mulby_coeff[0], lmax_norm2_mulby_coeff[1]);

    // * A = currentbestVectorLength - partial_sum_norm2 - (max_remaining_norm2 * remaining_coeff_2)
    bigring sub_of_abc;
    bigring bestVectorLength_big;
    ring_to_bigring(bestVectorLength_big, *bestVectorLength);
    bigring_sub(sub_of_abc, bestVectorLength_big, partial_sum_norm2_big);
    bigring_sub(sub_of_abc, sub_of_abc, lmax_norm2_mulby_coeff);

    // printf("sub_of_abc: (%lld, %lld)\n", sub_of_abc[0], sub_of_abc[1]);

    bigring zero = {0, 0};
    if (bigring_comp(sub_of_abc, zero) < 0) {
        return true;
    }

    // * B = 4 * (max_remaining_norm2 * remaining_coeff_2) * partial_sum_norm2
    bigring B;
    bigring_mul(B, lmax_norm2_mulby_coeff, partial_sum_norm2_big);
    bigring_scale(B, B, 4, 4);

    // printf("B: (%lld, %lld)\n", B[0], B[1]);

    // * C = A^2
    bigring sub_of_abc_squared;
    bigring_mul(sub_of_abc_squared, sub_of_abc, sub_of_abc);

    // printf("sub_of_abc_squared: (%lld, %lld)\n", sub_of_abc_squared[0], sub_of_abc_squared[1]);

    return bigring_comp(B, sub_of_abc_squared) > 0;
}

void update_best_available_index(bitmask_t used_mask, int numPoints) {
    while (best_available_index < numPoints / 2 && IS_USED(used_mask, best_available_index)) {
        ++best_available_index;
    }
}

static void find_best_permutation(int depth, const int numPoints, const vect *points, ring *bestVectorLength,
                                  int *bestPermutation, vect currentScaledVectorAtDepth, int *VArray,
                                  int *currentPermutation, bitmask_t used_mask, SymmetryMap *root,
                                  MyArray *centralSymmetryList, int32_t *coeff_sum_on_depths, const ring *pointNorms) {
    if (depth == numPoints / 2) {
        evaluate_permutation(currentScaledVectorAtDepth, numPoints, bestVectorLength, bestPermutation,
                             currentPermutation);
        return;
    }

    for (int candidate = 0; candidate < numPoints; ++candidate) {
        int mirror = centralSymmetryList->data[candidate];
        if (IS_USED(used_mask, candidate) || IS_USED(used_mask, mirror)) continue;

        SymmetryMap *nextChild = getChild(root, candidate);
        if (root != NULL && root->children_count > 0 && nextChild == NULL) continue;

        vect nextVectScaled;
        vect_scale(nextVectScaled, points[candidate], VArray[depth], currentScaledVectorAtDepth);
        if ((*bestVectorLength)[0] != -1 || (*bestVectorLength)[1] != -1) {
            if (!cheap_upper_bound_estimate(depth, nextVectScaled, bestVectorLength, coeff_sum_on_depths, pointNorms)) {
                continue;
            }

            // if (!should_check_branch_ub(depth, nextVectScaled, bestVectorLength, pointNorms, coeff_sum_on_depths)) {
            //     continue;  // Skip this branch if upper bound is not promising
            // }
        }

        currentPermutation[depth] = candidate;
        currentPermutation[numPoints - 1 - depth] = mirror;

        SET_USED(used_mask, candidate);
        SET_USED(used_mask, mirror);
        if (candidate == best_available_index) {
            update_best_available_index(used_mask, numPoints);
        }

        find_best_permutation(depth + 1, numPoints, points, bestVectorLength, bestPermutation, nextVectScaled, VArray,
                              currentPermutation, used_mask, nextChild, centralSymmetryList, coeff_sum_on_depths,
                              pointNorms);

        UNSET_USED(used_mask, candidate);
        UNSET_USED(used_mask, mirror);
        if (candidate < best_available_index) {
            best_available_index = candidate;
        }
    }
}

static void prepare_symmetry_graph(const int numPoints, const int *symIndexListToSearch,
                                   const int symIndexListToSearchLength, bool *alreadyUsedNodes, int *symmetryGroup,
                                   SymmetryMap *root) {
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

        SymmetryMap *child = createNode();

        addChild(root, nextChosen, child);

        prepare_symmetry_graph(numPoints, nextSymIndexList.data, nextSymIndexList.size, alreadyUsedNodes, symmetryGroup,
                               child);

        alreadyUsedNodes[nextChosen] = false;
    }
}

static void get_central_symmetry(const int numPoints, const vect *points, MyArray *out) {
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
    ring pointNorms[MAX_POINTS];
    int32_t coeff_sum_on_depths[MAX_POINTS] = {0};
    bitmask_t used_mask = 0;
    sqrt_our_k = sqrt((double)our_k);
    vect currentScaledVectorAtDepth = {{0, 0}, {0, 0}, {0, 0}};
    bool alreadyUsedNodes[MAX_POINTS] = {false};
    int currentPermutation[MAX_POINTS] = {0};
    int symmetryGroup[MAX_POINTS * MAX_SYM_GROUP_SIZE];
    precompute_point_norms(points, numPoints, pointNorms);

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

    printf("point norms precomputed\n");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", pointNorms[i][0]);
        printf("%d ", pointNorms[i][1]);
        printf("\n");
    }

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

    printf("coeff_sum_on_depths: ");
    for (int i = 0; i < numPoints / 2; ++i) {
        printf("%d ", coeff_sum_on_depths[i]);
    }
    printf("\n");

    MyArray centralSymmetryList;
    get_central_symmetry(numPoints, points, &centralSymmetryList);
    printf("Is centrally symmetric: %d\n", centralSymmetryList.size != 0);

    SymmetryMap *root = createNode();
    prepare_symmetry_graph(numPoints, symIndexListToSearch, symmetryGroupSize, alreadyUsedNodes, symmetryGroup, root);
    printf("symmetry map prep done\n");

    clock_t begin = clock();
    find_best_permutation(0, numPoints, points, bestVectorLength, bestPermutation, currentScaledVectorAtDepth, VArray,
                          currentPermutation, used_mask, root, &centralSymmetryList, coeff_sum_on_depths, pointNorms);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time spent: %f\n", time_spent);

    (*bestVectorLength)[0] = (*bestVectorLength)[0] * 4;
    (*bestVectorLength)[1] = (*bestVectorLength)[1] * 4;
}