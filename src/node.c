#include "node.h"
#include <math.h>

int find_or_add_node(Point p, Node *nodes, int *n_nodes) {
    for (int i = 0; i < *n_nodes; i++) {
        if (points_are_equal(nodes[i].position, p)) {
            return nodes[i].id;
        }
    }
    nodes[*n_nodes] = (Node){*n_nodes + 1, 0, p};
    return ++(*n_nodes);
}

void find_and_update(int *nA, int *nB, Point *body, int n_points, Node **nodes, int *n_nodes) {
    Point p;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;
        if (get_intersection((*nodes)[*nA-1].position, (*nodes)[*nB-1].position, body[i], body[j], &p)) {
            break;
        }
    }
    *nA = find_or_add_node(p, *nodes, n_nodes);
}