#ifndef ELEMENT_H
#define ELEMENT_H

#include "node.h"

typedef struct {
    int id;
    int type;  // 2 = triangle, 3 = quadrilateral, 4 = pentagon
    int n_nodes; // 3 = triangle, 4 = quadrilateral, 5 = pentagon
    int nodes[5];
    int flag; // n1=1 n2=2 n3=3 n4=4
} Element;

void tria_vertices(int *n1, int *n2, int *n3, int *n4, Point *body, int n_points, Node **nodes, int *n_nodes);
void quad_vertices(int *n1, int *n2, int *n3, int *n4, Point *body, int n_points, Node **nodes, int *n_nodes);
int penta_vertices(int *n1, int *n2, int *n3, int *n4, int *n5, Point *body, int n_points, Node **nodes, int *n_nodes);
void split_pentagons(Element **elements, int *num_elements);

#endif