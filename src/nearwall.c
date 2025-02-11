#include "nearwall.h"
#include "point.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void compute_offset(Point **offset, Point *body, int n_points, double nwl_distance) {

    *offset = (Point *)malloc(n_points*sizeof(Point));

    double A = 0.0;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;
        A += (body[i].x * body[j].y) - (body[j].x * body[i].y);
    }
    int v = 2 * (A < 0) - 1;

    for (int i = 0; i < n_points; i++) {
        int h = (n_points + i - 1) % n_points;
        int j = (i + 1) % n_points;

        double hj = sqrt((body[j].x-body[h].x)*(body[j].x-body[h].x) + (body[j].y-body[h].y)*(body[j].y-body[h].y));
        double nx = - v * (body[j].y - body[h].y)/hj;
        double ny = + v * (body[j].x - body[h].x)/hj;

        (*offset)[i].x = body[i].x + nx * nwl_distance;
        (*offset)[i].y = body[i].y + ny * nwl_distance;
    }

    // TODO: OFFSET QUALITY CHECKS
}

double get_distance(Node a, Node b) {
    double dx = a.position.x - b.position.x;
    double dy = a.position.y - b.position.y;
    return sqrt(dx*dx + dy*dy);
}

int find_nearest(Node *nodes, int *ids, int n, int current, int *flag, double cell_size) {
    int nearest = -1;
    double min_distance = 2 * cell_size;

    for (int i = 0; i < n; i++) {
        if (!flag[i]) {
            double distance = get_distance(nodes[current], nodes[i]);
            if (distance < min_distance) {
                min_distance = distance;
                nearest = i;
            }
        }
    }
    return nearest;
}

void get_offset_nodes(int **offset_nodes, int *n_offset_nodes, Node *nodes, int n_nodes, double cell_size) {

    for (int i = 0; i < n_nodes; i++) {
        if (nodes[i].type == 2) (*n_offset_nodes)++;
    }

    *offset_nodes = malloc(*n_offset_nodes * sizeof(int));
    if (!offset_nodes) {
        perror("An error has occurred");
        return;
    }

    Node type2_nodes[*n_offset_nodes];
    int ids[*n_offset_nodes];
    int count = 0;
    for (int i = 0; i < n_nodes; i++) {
        if (nodes[i].type == 2) {
            type2_nodes[count] = nodes[i];
            ids[count] = nodes[i].id;
            count++;
        }
    }

    int flag[*n_offset_nodes];
    for (int i = 0; i < *n_offset_nodes; i++) flag[i] = 0;

    int current = 0;
    flag[current] = 1;
    (*offset_nodes)[0] = ids[current];

    for (int i = 1; i < *n_offset_nodes; i++) {
        int nearest = find_nearest(type2_nodes, ids, *n_offset_nodes, current, flag, cell_size);
        if (nearest == -1) break;

        flag[nearest] = 1;
        (*offset_nodes)[i] = ids[nearest];
        current = nearest;
    }
}

void extrude_near_wall_cells(
    Element **elements, int *n_elements,
    Node **nodes, int *n_nodes,
    Point *body, int n_points,
    int *offset_nodes, int n_offset_nodes,
    double nwl_distance
) {
    Point pA[n_offset_nodes];
    for (int i = 0; i < n_offset_nodes; i++) {
        pA[i] = (*nodes)[offset_nodes[i]-1].position;
    }

    double A = 0.0;
    for (int i = 0; i < n_offset_nodes; i++) {
        int j = (i + 1) % n_offset_nodes;
        A += (pA[i].x * pA[j].y) - (pA[j].x * pA[i].y);
    }
    int v = 2 * (A > 0) - 1;

    Point pB[n_offset_nodes];
    for (int i = 0; i < n_offset_nodes; i++) {
        int h = (n_offset_nodes + i - 1) % n_offset_nodes;
        int j = (i + 1) % n_offset_nodes;

        double hj = sqrt((pA[j].x-pA[h].x)*(pA[j].x-pA[h].x) + (pA[j].y-pA[h].y)*(pA[j].y-pA[h].y));
        double nx = - v * (pA[j].y - pA[h].y)/hj;
        double ny = + v * (pA[j].x - pA[h].x)/hj;

        pB[i].x = pA[i].x + nx * nwl_distance;
        pB[i].y = pA[i].y + ny * nwl_distance;
    }

    Point pI[n_offset_nodes];
    for (int k = 0; k < n_offset_nodes; k++) {
        c = 0;
        for (int i = 0; i < n_points; i++) {
            int j = (i + 1) % n_points;
            if (get_intersection(pA[k], pB[k], body[i], body[j], &(pI[k]))) {
                break;
            }
        }
    }

    for (int i = 0; i < n_offset_nodes; i++) {
        printf("%d ", offset_nodes[i]);
    }

    // CREATE NODES ON INTERSECTION POINTS
    // CREATE NODES BASED ON DISTRIBUTION
    // CREATE ELEMENTS
}