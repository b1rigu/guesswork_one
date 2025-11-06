#include <stdbool.h>
#include <stdckdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ring128.h"
#include "symmetry.h"
#include "utils.h"

typedef struct {
    ring norm;
    int original_index;
} NormRingWithIndex;

NormRingWithIndex sortedRings[MAX_POINTS];

int compare_ring_desc(const void* a, const void* b) {
    const NormRingWithIndex* ra = (const NormRingWithIndex*)a;
    const NormRingWithIndex* rb = (const NormRingWithIndex*)b;
    return -ring_comp(ra->norm, rb->norm);  // descending order
}

void precompute_sorted_norms(const vect* points, int numPoints, NormRingWithIndex* outRings) {
    for (int i = 0; i < numPoints; ++i) {
        vect_norm2(outRings[i].norm, points[i]);
        outRings[i].original_index = i;
    }

    qsort(outRings, numPoints, sizeof(NormRingWithIndex), compare_ring_desc);
}

const ring* get_best_unused_ring(uint64_t used_mask, const NormRingWithIndex* sortedRings, int numPoints) {
    for (int i = 0; i < numPoints; ++i) {
        int idx = sortedRings[i].original_index;
        if (!IS_USED(used_mask, idx)) {
            return &sortedRings[i].norm;
        }
    }
    return NULL;
}

static bool should_check_branch_ub(int depth, const vect nextScaledVect, const ring* bestVectorLength,
                                   int32_t* coeff_sum2_on_depths, bitmask_t used_mask, int numPoints) {
    // ! A = currentbestVectorLength - partial_sum_norm2 - (max_remaining_norm2 * remaining_coeff_2)
    // ! B = 4 * (max_remaining_norm2 * remaining_coeff_2) * partial_sum_norm2
    // ! C = A^2
    ring partial_sum_norm2;  // * a^2
    vect_norm2(partial_sum_norm2, nextScaledVect);
    bigring partial_sum_norm2_big;
    ring_to_bigring(partial_sum_norm2_big, partial_sum_norm2);

    int64_t remaining_coeff_2 = coeff_sum2_on_depths[depth + 1];

    const ring* current_max_norm2 = get_best_unused_ring(used_mask, sortedRings, numPoints);

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

static void evaluate_permutation(vect currentScaledVector, int numPoints, ring* bestVectorLength, int* bestPermutation,
                                 int* currentPermutation) {
    ring currentSquared;
    vect_norm2(currentSquared, currentScaledVector);
    if (ring_comp(currentSquared, *bestVectorLength) > 0) {
        (*bestVectorLength)[0] = currentSquared[0];
        (*bestVectorLength)[1] = currentSquared[1];
        memcpy(bestPermutation, currentPermutation, numPoints * sizeof(*currentPermutation));
    }
}

// * centralSymmetryList->size != 0 equals to the states having central symmetry
static void find_best_permutation(int depth, const int numPoints, const vect* points, ring* bestVectorLength,
                                  int* bestPermutation, vect currentScaledVectorAtDepth, int* VArray,
                                  int* currentPermutation, bitmask_t used_mask, SymmetryMap* root,
                                  MyArray* centralSymmetryList, int32_t* coeff_sum2_on_depths) {
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

        SymmetryMap* nextChild = getChild(root, candidate);
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

static void get_central_symmetry(const int numPoints, const vect* points, MyArray* out) {
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

static void generate_zero_sum_array(int numPoints, int* VArray) {
    int start = numPoints - 1;
    for (int i = 0; i < numPoints; ++i) {
        VArray[i] = start - 2 * i;
    }
}

int main() {
    vect points[MAX_POINTS];
    int numPoints = nextInt();
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < 3; j++) {
            points[i][j][0] = nextLL();  // a part
            points[i][j][1] = nextLL();  // b part
        }
    }

    int bestPermutation[MAX_POINTS];
    int VArray[MAX_POINTS];

    generate_zero_sum_array(numPoints, VArray);
    precompute_sorted_norms(points, numPoints, sortedRings);

    int32_t coeff_sum2_on_depths[MAX_POINTS] = {0};
    int total = 0;
    for (int i = numPoints / 2 - 1; i >= 0; --i) {
        total += VArray[i];
        coeff_sum2_on_depths[i] = total * total;
    }
    for (int i = 0; i < numPoints / 2; ++i) {
        coeff_sum2_on_depths[numPoints - 1 - i] = coeff_sum2_on_depths[i];
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

    bitmask_t used_mask = 0;
    vect currentScaledVectorAtDepth = {{0, 0}, {0, 0}, {0, 0}};
    int currentPermutation[MAX_POINTS] = {0};

    MyArray centralSymmetryList;
    get_central_symmetry(numPoints, points, &centralSymmetryList);
    printf("Is centrally symmetric: %d\n", centralSymmetryList.size != 0);

    ring bestVectorLength = {-1, -1};

    SymmetryMap* root = get_symmetry_explore_map(numPoints, points);

    clock_t begin = clock();
    find_best_permutation(0, numPoints, points, &bestVectorLength, bestPermutation, currentScaledVectorAtDepth, VArray,
                          currentPermutation, used_mask, root, &centralSymmetryList, coeff_sum2_on_depths);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time spent: %f\n", time_spent);

    bestVectorLength[0] = bestVectorLength[0] * 4;
    bestVectorLength[1] = bestVectorLength[1] * 4;

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