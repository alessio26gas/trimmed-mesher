#include "nearwall.h"
#include <math.h>
#include <stdlib.h>

void compute_offset(Point **offset, Point *body, int n_points, double nwl_distance) {

    *offset = (Point *)malloc(n_points*sizeof(Point));

    for (int i = 0; i < n_points; i++) {
        int h = (n_points + i - 1) % n_points;
        int j = (i + 1) % n_points;

        double hj = sqrt((body[j].x-body[h].x)*(body[j].x-body[h].x) + (body[j].y-body[h].y)*(body[j].y-body[h].y));

        int v = -1; // TODO

        double nx = - v * (body[j].y - body[h].y)/hj;
        double ny = + v * (body[j].x - body[h].x)/hj;

        (*offset)[i].x = body[i].x + nx * nwl_distance;
        (*offset)[i].y = body[i].y + ny * nwl_distance;
    }
}