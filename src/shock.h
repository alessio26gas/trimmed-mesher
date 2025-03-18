#ifndef SHOCK_H
#define SHOCK_H

#include "nearwall.h"

typedef struct {
    double distance;
    double first;
    double last;
    double SF;
    int n;
    int distribution;
} NearShockLayer;

void compute_nsl(NearShockLayer *nsl);
bool compute_shock_offset(Point **offset, int *n_offset, Point *shock, int n_shock, double nsl_distance, double cell_size);
void extrude_near_shock_cells(
    Element **elements, int *n_elements,
    Node **nodes, int *n_nodes,
    Point *shock, int n_shock,
    int *offset_nodes, int n_offset_nodes,
    NearShockLayer nsl, int simmetry,
    double X0, int cols, double cell_size
);

#endif