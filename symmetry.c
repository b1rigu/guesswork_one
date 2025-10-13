#include "symmetry.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ring128.h"
#include "utils.h"

typedef struct {
    int data[MAX_SYM_GROUP_SIZE];
    int size;
} MyArray;

typedef struct SymmetryMap {
    signed char children_index_list[MAX_POINTS];
    unsigned char children_count;
    struct SymmetryMap **children;
} SymmetryMap;

typedef struct {
    ring norm;
    int original_index;
} NormRingWithIndex;

NormRingWithIndex sortedRings[MAX_POINTS];

typedef __uint64_t bitmask_t;

#define IS_USED(mask, i) (((mask) >> (i)) & 1)
#define SET_USED(mask, i) ((mask) |= ((bitmask_t)1 << (i)))
#define UNSET_USED(mask, i) ((mask) &= ~((bitmask_t)1 << (i)))

// row 0 :  table[0*N + 0]  table[0*N + 1] … table[0*N + N-1]
// row 1 :  table[1*N + 0]  table[1*N + 1] … table[1*N + N-1]
// row 2 :  ...
#define SYM_IMAGE(table, row, N, column) ((table)[(size_t)(row) * (N) + (column)])

static bool is_linear_dependant(const int n, const vect point1, const vect point2, const vect point3) {
    if (n <= 1) {
        return false;
    } else if (n == 2) {
        ring sum, temp1, temp2, temp3;
        ring_mul(temp1, point1[0], point1[0]);
        ring_mul(temp2, point1[1], point1[1]);
        ring_mul(temp3, point1[2], point1[2]);
        ring_add(sum, temp1, temp2);
        ring_add(sum, sum, temp3);

        ring sum2;
        ring_mul(temp1, point2[0], point2[0]);
        ring_mul(temp2, point2[1], point2[1]);
        ring_mul(temp3, point2[2], point2[2]);
        ring_add(sum2, temp1, temp2);
        ring_add(sum2, sum2, temp3);

        ring left;
        ring_mul(left, sum, sum2);

        ring_mul(temp1, point1[0], point2[0]);
        ring_mul(temp2, point1[1], point2[1]);
        ring_mul(temp3, point1[2], point2[2]);
        ring_add(sum, temp1, temp2);
        ring_add(sum, sum, temp3);

        ring right;
        ring_mul(right, sum, sum);

        ring total;
        ring_sub(total, left, right);

        ring zero = {0, 0};
        if (ring_comp(total, zero) == 0) {
            return true;
        }
        return false;
    } else if (n == 3) {
        ring first, second, third, fourth, fifth, sixth, temp1, temp2;

        ring_mul(temp1, point1[0], point2[1]);
        ring_mul(first, temp1, point3[2]);

        ring_mul(temp1, point1[1], point2[2]);
        ring_mul(second, temp1, point3[0]);

        ring_mul(temp1, point2[0], point3[1]);
        ring_mul(third, temp1, point1[2]);

        ring_mul(temp1, point1[2], point2[1]);
        ring_mul(fourth, temp1, point3[0]);

        ring_mul(temp1, point1[1], point2[0]);
        ring_mul(fifth, temp1, point3[2]);

        ring_mul(temp1, point2[2], point3[1]);
        ring_mul(sixth, temp1, point1[0]);

        ring total;
        ring_add(total, first, second);
        ring_add(total, total, third);
        ring_sub(total, total, fourth);
        ring_sub(total, total, fifth);
        ring_sub(total, total, sixth);

        ring zero = {0, 0};
        if (ring_comp(total, zero) == 0) {
            return true;
        }
        return false;
    }

    return true;
}

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
        printf("Tried to overwrite child at pos: %d\n", pos);
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

int compare_ring_desc(const void *a, const void *b) {
    const NormRingWithIndex *ra = (const NormRingWithIndex *)a;
    const NormRingWithIndex *rb = (const NormRingWithIndex *)b;
    return -ring_comp(ra->norm, rb->norm);  // descending order
}

void precompute_sorted_norms(const vect *points, int numPoints, NormRingWithIndex *outRings) {
    for (int i = 0; i < numPoints; ++i) {
        vect_norm2(outRings[i].norm, points[i]);
        outRings[i].original_index = i;
    }

    qsort(outRings, numPoints, sizeof(NormRingWithIndex), compare_ring_desc);
}

const ring *get_best_unused_ring(uint64_t used_mask, const NormRingWithIndex *sortedRings, int numPoints) {
    for (int i = 0; i < numPoints; ++i) {
        int idx = sortedRings[i].original_index;
        if (!IS_USED(used_mask, idx)) {
            return &sortedRings[i].norm;
        }
    }
    return NULL;
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

static bool should_check_branch_ub(int depth, const vect nextScaledVect, const ring *bestVectorLength,
                                   int32_t *coeff_sum2_on_depths, bitmask_t used_mask, int numPoints) {
    // ! A = currentbestVectorLength - partial_sum_norm2 - (max_remaining_norm2 * remaining_coeff_2)
    // ! B = 4 * (max_remaining_norm2 * remaining_coeff_2) * partial_sum_norm2
    // ! C = A^2
    ring partial_sum_norm2;  // * a^2
    vect_norm2(partial_sum_norm2, nextScaledVect);
    bigring partial_sum_norm2_big;
    ring_to_bigring(partial_sum_norm2_big, partial_sum_norm2);

    int64_t remaining_coeff_2 = coeff_sum2_on_depths[depth + 1];

    const ring *current_max_norm2 = get_best_unused_ring(used_mask, sortedRings, numPoints);

    bigring lmax_norm2_mulby_coeff;
    bigring current_max_norm2_big;
    ring_to_bigring(current_max_norm2_big, *current_max_norm2);
    bigring_scale(lmax_norm2_mulby_coeff, current_max_norm2_big, remaining_coeff_2, remaining_coeff_2);

    // * A = currentbestVectorLength - partial_sum_norm2 - (max_remaining_norm2 * remaining_coeff_2)
    bigring sub_of_abc;
    bigring bestVectorLength_big;
    ring_to_bigring(bestVectorLength_big, *bestVectorLength);
    bigring_sub(sub_of_abc, bestVectorLength_big, partial_sum_norm2_big);
    bigring_sub(sub_of_abc, sub_of_abc, lmax_norm2_mulby_coeff);

    bigring zero = {0, 0};
    if (bigring_comp(sub_of_abc, zero) < 0) {
        return true;
    }

    // * B = 4 * (max_remaining_norm2 * remaining_coeff_2) * partial_sum_norm2
    bigring B;
    bigring_mul(B, lmax_norm2_mulby_coeff, partial_sum_norm2_big);
    bigring_scale(B, B, 4, 4);

    // * C = A^2
    bigring sub_of_abc_squared;
    bigring_mul(sub_of_abc_squared, sub_of_abc, sub_of_abc);

    return bigring_comp(B, sub_of_abc_squared) > 0;
}

static void find_best_permutation(int depth, const int numPoints, const vect *points, ring *bestVectorLength,
                                  int *bestPermutation, vect currentScaledVectorAtDepth, int *VArray,
                                  int *currentPermutation, bitmask_t used_mask, SymmetryMap *root,
                                  MyArray *centralSymmetryList, int32_t *coeff_sum2_on_depths) {
    if (depth == numPoints || (centralSymmetryList->size != 0 && depth == numPoints / 2)) {
        evaluate_permutation(currentScaledVectorAtDepth, numPoints, bestVectorLength, bestPermutation,
                             currentPermutation);
        return;
    }

    for (int candidate = 0; candidate < numPoints; ++candidate) {
        if (IS_USED(used_mask, candidate)) continue;
        int mirror = -1;
        if (centralSymmetryList->size != 0) {
            mirror = centralSymmetryList->data[candidate];
            if (IS_USED(used_mask, mirror)) continue;
        }

        SymmetryMap *nextChild = getChild(root, candidate);
        if (root != NULL && root->children_count > 0 && nextChild == NULL) continue;

        vect nextVectScaled;
        vect_scale(nextVectScaled, points[candidate], VArray[depth], currentScaledVectorAtDepth);
        if ((*bestVectorLength)[0] != -1 || (*bestVectorLength)[1] != -1) {
            if (!should_check_branch_ub(depth, nextVectScaled, bestVectorLength, coeff_sum2_on_depths, used_mask,
                                        numPoints)) {
                continue;
            }
        }

        currentPermutation[depth] = candidate;
        if (mirror != -1) {
            currentPermutation[numPoints - 1 - depth] = mirror;
        }

        SET_USED(used_mask, candidate);
        if (mirror != -1) {
            SET_USED(used_mask, mirror);
        }

        find_best_permutation(depth + 1, numPoints, points, bestVectorLength, bestPermutation, nextVectScaled, VArray,
                              currentPermutation, used_mask, nextChild, centralSymmetryList, coeff_sum2_on_depths);

        UNSET_USED(used_mask, candidate);
        if (mirror != -1) {
            UNSET_USED(used_mask, mirror);
        }
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

static void prepare_symmetry_graph(const int numPoints, const int *symIndexListToSearch,
                                   const int symIndexListToSearchLength, bool *alreadyUsedNodes, int *symmetryGroup,
                                   SymmetryMap *root, int depth, const vect *points, int *permutationTracker) {
    MyArray nodes_to_explore;
    get_next_nodes_to_explore(numPoints, symIndexListToSearch, symIndexListToSearchLength, &nodes_to_explore,
                              alreadyUsedNodes, symmetryGroup);

    for (int i = 0; i < nodes_to_explore.size; ++i) {
        int nextChosen = nodes_to_explore.data[i];
        alreadyUsedNodes[nextChosen] = true;

        if (depth == 1 && is_linear_dependant(2, points[permutationTracker[0]], points[nextChosen], NULL)) {
            addChild(root, nextChosen, root);
        } else if (depth == 2 && is_linear_dependant(3, points[permutationTracker[0]], points[permutationTracker[1]],
                                                     points[nextChosen])) {
            addChild(root, nextChosen, root);
        } else if (depth == 3) {
            addChild(root, nextChosen, root);
        } else {
            MyArray nextSymIndexList;
            get_next_symmetry_group(symIndexListToSearch, symIndexListToSearchLength, numPoints, nextChosen,
                                    &nextSymIndexList, symmetryGroup);

            permutationTracker[depth] = nextChosen;

            SymmetryMap *child = createNode();

            addChild(root, nextChosen, child);

            prepare_symmetry_graph(numPoints, nextSymIndexList.data, nextSymIndexList.size, alreadyUsedNodes,
                                   symmetryGroup, child, depth + 1, points, permutationTracker);
        }

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

static void generate_zero_sum_array(int numPoints, int *VArray) {
    int start = (numPoints - 1) * 1;  // largest positive value
    for (int i = 0; i < numPoints; ++i) {
        VArray[i] = start - 2 * i;
    }
}

void use_symmetry(int numPoints, vect *points, ring *bestVectorLength, int *bestPermutation) {
    int symmetryGroup[MAX_POINTS * MAX_SYM_GROUP_SIZE];
    int symmetryGroupSize = nextInt();
    for (int g = 0; g < symmetryGroupSize; g++) {
        for (int i = 0; i < numPoints; i++) {
            symmetryGroup[g * numPoints + i] = nextInt();
        }
    }

    int VArray[MAX_POINTS];
    generate_zero_sum_array(numPoints, VArray);
    
    precompute_sorted_norms(points, numPoints, sortedRings);

    int32_t coeff_sum2_on_depths[MAX_POINTS] = {0};
    int total = 0;
    for (int i = numPoints / 2 - 1; i >= 0; --i) {
        total += VArray[i];
        coeff_sum2_on_depths[i] = total * total;
    }

    printf("varray: ");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", VArray[i]);
    }
    printf("\n");

    printf("coeff_sum2_on_depths: ");
    for (int i = 0; i < numPoints; ++i) {
        printf("%d ", coeff_sum2_on_depths[i]);
    }
    printf("\n");

    int symIndexListToSearch[MAX_SYM_GROUP_SIZE];
    for (int g = 0; g < symmetryGroupSize; ++g) {
        symIndexListToSearch[g] = g;
    }

    bitmask_t used_mask = 0;
    vect currentScaledVectorAtDepth = {{0, 0}, {0, 0}, {0, 0}};
    bool alreadyUsedNodes[MAX_POINTS] = {false};
    int currentPermutation[MAX_POINTS] = {0};

    MyArray centralSymmetryList;
    get_central_symmetry(numPoints, points, &centralSymmetryList);
    printf("Is centrally symmetric: %d\n", centralSymmetryList.size != 0);

    int permutationTracker[MAX_POINTS] = {-1};
    SymmetryMap *root = createNode();
    prepare_symmetry_graph(numPoints, symIndexListToSearch, symmetryGroupSize, alreadyUsedNodes, symmetryGroup, root, 0,
                           points, permutationTracker);
    printf("symmetry map prep done\n");

    clock_t begin = clock();
    find_best_permutation(0, numPoints, points, bestVectorLength, bestPermutation, currentScaledVectorAtDepth, VArray,
                          currentPermutation, used_mask, root, &centralSymmetryList, coeff_sum2_on_depths);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time spent: %f\n", time_spent);

    (*bestVectorLength)[0] = (*bestVectorLength)[0] * 4;
    (*bestVectorLength)[1] = (*bestVectorLength)[1] * 4;
}