#ifndef NEARWALL_H
#define NEARWALL_H

#include "element.h"

typedef struct {
    double distance;
    double first;
    double last;
    double SF;
    int n;
    double min_surf_distance;
    int surf_max_iter;
    int distribution;
} NearWallLayer;

int get_nwl_n(NearWallLayer nwl);
double get_nwl_distance(NearWallLayer nwl);
double get_SF(NearWallLayer nwl);
void compute_nwl(NearWallLayer *nwl);
bool compute_offset(Point **offset, int *n_offset, Point *body, int n_body, double nwl_distance, double cell_size);
int get_offset_nodes(int **offset_nodes, int *n_offset_nodes, Node *nodes, int n_nodes, double cell_size, double X0, double Y0, int rows, int cols);
void extrude_near_wall_cells(
    Element **elements, int *n_elements,
    Node **nodes, int *n_nodes,
    Point *body, int n_body,
    int *offset_nodes, int n_offset_nodes,
    NearWallLayer nwl, int simmetry,
    double X0, int cols, double cell_size
);

#endif