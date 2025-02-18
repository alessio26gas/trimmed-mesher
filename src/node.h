#ifndef NODE_H
#define NODE_H

#include "point.h"

typedef struct {
    int id;
    int type; // 0 = not enclosed, 1 = enclosed, 2 = near-body
    Point position;
} Node;

void find_and_update(int *nA, int *nB, Point *body, int n_points, Node **nodes, int *n_nodes);

#endif