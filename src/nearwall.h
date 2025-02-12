#ifndef NEARWALL_H
#define NEARWALL_H

#include "element.h"

typedef struct {
    double distance;
    double first;
    double last;
    double SF;
    int n;
} NearWallLayer;

double get_SF(NearWallLayer nwl);
bool compute_offset(Point **offset, int *n_offset, Point *body, int n_body, double nwl_distance);
void get_offset_nodes(int **offset_nodes, int *n_offset_nodes, Node *nodes, int n_nodes, double cell_size);
void extrude_near_wall_cells(
    Element **elements, int *n_elements,
    Node **nodes, int *n_nodes,
    Point *body, int n_body,
    int *offset_nodes, int n_offset_nodes,
    NearWallLayer nwl
);

#endif