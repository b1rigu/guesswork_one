#include "symmetry.h"

#include <math.h>
#include <stdbool.h>
#include <stdckdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct SymmetryMap {
    int indexListToSearch[MAX_POINTS];
    int indexListSize;
    struct SymmetryMap *children[MAX_POINTS];
} SymmetryMap;

SymmetryMap allNodes[MAX_POINTS * MAX_POINTS];
int nextFreeNode = 0;

int program_subloops = 0;
int program_mainloops = 0;

// row 0 :  table[0*N + 0]  table[0*N + 1] … table[0*N + N-1]
// row 1 :  table[1*N + 0]  table[1*N + 1] … table[1*N + N-1]
// row 2 :  ...
#define SYM_IMAGE(table, row, N, column) ((table)[(size_t)(row) * (N) + (column)])

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

static void find_best_permutation(int depth, const int numPoints, const vect *points, ring *bestVectorLength,
                                  int *bestPermutation, vect currentScaledVectorAtDepth, int *VArray,
                                  int *currentPermutation, bool *alreadyUsedNodes, SymmetryMap *root) {
    program_mainloops++;
    if (depth == numPoints) {
        evaluate_permutation(currentScaledVectorAtDepth, numPoints, bestVectorLength, bestPermutation,
                             currentPermutation);
        return;
    }

    for (int candidate = 0; candidate < numPoints; ++candidate) {
        if (alreadyUsedNodes[candidate]) continue;
        if (root != NULL && root->indexListSize > 0 && root->children[candidate] == NULL) continue;

        alreadyUsedNodes[candidate] = true;
        currentPermutation[depth] = candidate;

        vect nextVectScaled;
        vect_scale(&nextVectScaled, points[candidate], VArray[depth], currentScaledVectorAtDepth);

        SymmetryMap *nextChild = (root != NULL && root->indexListSize > 0) ? root->children[candidate] : NULL;

        find_best_permutation(depth + 1, numPoints, points, bestVectorLength, bestPermutation, nextVectScaled, VArray,
                              currentPermutation, alreadyUsedNodes, nextChild);

        alreadyUsedNodes[candidate] = false;
    }
}

void initSymmetryMap(SymmetryMap *map) {
    map->indexListSize = 0;
    for (int i = 0; i < MAX_POINTS; i++) {
        map->children[i] = NULL;
    }
}

SymmetryMap *getNewNode() {
    if (nextFreeNode >= MAX_POINTS * MAX_POINTS) {
        // Error: out of nodes
        return NULL;
    }

    SymmetryMap *newNode = &allNodes[nextFreeNode++];
    initSymmetryMap(newNode);
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
        root->children[nextChosen] = child;
        root->indexListToSearch[root->indexListSize++] = nextChosen;

        prepare_symmetry_graph(numPoints, nextSymIndexList.data, nextSymIndexList.size, alreadyUsedNodes, symmetryGroup,
                               child);

        alreadyUsedNodes[nextChosen] = false;
    }
}

void use_symmetry(int numPoints, vect *points, ring *bestVectorLength, int *bestPermutation, int *VArray) {
    bool alreadyUsedNodes[MAX_POINTS] = {false};
    int currentPermutation[MAX_POINTS] = {0};
    int symmetryGroup[MAX_POINTS * MAX_SYM_GROUP_SIZE];

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

    vect currentScaledVectorAtDepth = {{0, 0}, {0, 0}, {0, 0}};

    SymmetryMap *root = getNewNode();
    prepare_symmetry_graph(numPoints, symIndexListToSearch, symmetryGroupSize, alreadyUsedNodes, symmetryGroup, root);

    memset(alreadyUsedNodes, false, sizeof(alreadyUsedNodes));

    find_best_permutation(0, numPoints, points, bestVectorLength, bestPermutation, currentScaledVectorAtDepth, VArray,
                          currentPermutation, alreadyUsedNodes, root);

    printf("Total Runtime Count subloops: %d\n", program_subloops);
    printf("Total Runtime Count mainloops: %d\n", program_mainloops);
}