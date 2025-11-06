#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "ring.h"
#include "vect.h"

#define MAX_POINTS 50
#define MAX_SYM_GROUP_SIZE 200

typedef struct SymmetryMap {
    signed char children_index_list[MAX_POINTS];
    unsigned char children_count;
    struct SymmetryMap **children;
} SymmetryMap;

typedef struct {
    int data[MAX_SYM_GROUP_SIZE];
    int size;
} MyArray;

SymmetryMap *getChild(SymmetryMap *node, unsigned char index);

SymmetryMap *get_symmetry_explore_map(int numPoints, vect *points);

#endif