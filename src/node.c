#include "node.h"
#include <math.h>

void find_and_update(int *nA, int *nB, Point *body, int n_points, Node **nodes, int *n_nodes) {
    Point p;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;
        if (get_intersection((*nodes)[*nA-1].position, (*nodes)[*nB-1].position, body[i], body[j], &p)) {
            break;
        }
    }
    (*nodes)[*nA - 1].position.x = p.x;
    (*nodes)[*nA - 1].position.y = p.y;
    (*nodes)[*nA - 1].type = 2;
}