#ifndef NEARWALL_H
#define NEARWALL_H

#include "nearwall.h"
#include "point.h"
#include "node.h"
#include "element.h"

void compute_offset(Point **offset, Point *body, int n_points, double nwl_distance);
void get_offset_nodes(int **offset_nodes, int *n_offset_nodes, Node *nodes, int n_nodes, double cell_size);
void extrude_near_wall_cells(
    Element **elements, int *n_elements,
    Node **nodes, int *n_nodes,
    Point *body, int n_points,
    int *offset_nodes, int n_offset_nodes,
    double nwl_distance
);

#endif