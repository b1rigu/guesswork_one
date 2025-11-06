#include "symmetry.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "ring128.h"
#include "utils.h"

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
        ring first, second, third, fourth, fifth, sixth, temp;

        ring_mul(temp, point1[0], point2[1]);
        ring_mul(first, temp, point3[2]);

        ring_mul(temp, point1[1], point2[2]);
        ring_mul(second, temp, point3[0]);

        ring_mul(temp, point2[0], point3[1]);
        ring_mul(third, temp, point1[2]);

        ring_mul(temp, point1[2], point2[1]);
        ring_mul(fourth, temp, point3[0]);

        ring_mul(temp, point1[1], point2[0]);
        ring_mul(fifth, temp, point3[2]);

        ring_mul(temp, point2[2], point3[1]);
        ring_mul(sixth, temp, point1[0]);

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

static SymmetryMap *createNode() {
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

static void addChild(SymmetryMap *node, unsigned char index, SymmetryMap *child) {
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

SymmetryMap *getChild(SymmetryMap *node, unsigned char index) {
    if (node == NULL) return NULL;
    signed char pos = node->children_index_list[index];
    if (pos == -1) return NULL;
    return node->children[pos];
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

SymmetryMap *get_symmetry_explore_map(int numPoints, vect *points) {
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
    
    bool alreadyUsedNodes[MAX_POINTS] = {false};
    int permutationTracker[MAX_POINTS] = {-1};
    SymmetryMap *root = createNode();
    prepare_symmetry_graph(numPoints, symIndexListToSearch, symmetryGroupSize, alreadyUsedNodes, symmetryGroup, root, 0,
                           points, permutationTracker);
    printf("symmetry map prep done\n");

    return root;
}