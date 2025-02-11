#ifndef POINT_H
#define POINT_H

#include <stdbool.h>

typedef struct {
    double x, y;
} Point;

int load_body(const char *filename, Point **body, int *n_points);
int get_point_type(Point p, Point *body, int n_points);
bool is_near_body(Point *p, Point *body, int n_points, double cell_size, double frac);
bool get_intersection(Point a, Point b, Point c, Point d, Point *intersect);
bool points_are_equal(Point a, Point b);

#endif