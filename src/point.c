#include "point.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define EPSILON 1e-8

int load_body(const char *filename, Point **body, int *n_points) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("An error has occured while opening the file");
        return -1;
    }

    int count = 0;
    double x, y;
    while (fscanf(file, "%lf,%lf", &x, &y) == 2) {
        count++;
    }

    *body = (Point *)malloc(count * sizeof(Point));
    if (!*body) {
        perror("An error has occurred");
        fclose(file);
        return -1;
    }

    rewind(file);
    int i = 0;
    while (fscanf(file, "%lf,%lf", &(*body)[i].x, &(*body)[i].y) == 2) {
        i++;
    }

    *n_points = count;
    fclose(file);
    return 0;
}

int is_enclosed(Point p, Point *body, int n_points) {
    int inside = 0;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;
        if (((body[i].y > p.y) != (body[j].y > p.y)) &&
            (p.x < (body[j].x - body[i].x) * (p.y - body[i].y) / (body[j].y - body[i].y) + body[i].x)) {
            inside = !inside;
        }
    }
    return inside;
}

bool is_near_body(Point *p, Point *body, int n_points, double cell_size, double frac) {
    double min_distance = cell_size;
    double ABx, ABy, APx, APy, ABAB, ABAP, t, Xq, Yq, dx, dy, distance;
    double xp = (*p).x, yp = (*p).y;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;

        ABx = body[j].x - body[i].x;
        ABy = body[j].y - body[i].y;
        APx = (*p).x - body[i].x;
        APy = (*p).y - body[i].y;
        t = (APx * ABx + APy * ABy)/(ABx * ABx + ABy * ABy);
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        Xq = body[i].x + t * ABx;
        Yq = body[i].y + t * ABy;
        dx = Xq - (*p).x;
        dy = Yq - (*p).y;
        distance = sqrt(dx * dx + dy * dy);

        if (distance < min_distance) {
            min_distance = distance;
            xp = Xq;
            yp = Yq;
        }
    }
    if (min_distance < cell_size / frac) {
        (*p).x = xp;
        (*p).y = yp;
        return true;
    }
    return false;
}

bool get_intersection(Point a, Point b, Point c, Point d, Point *intersect) {
    double det = (b.x - a.x) * (d.y - c.y) - (b.y - a.y) * (d.x - c.x);
    if (fabs(det) < EPSILON) return false;
    double t = ((c.x - a.x) * (d.y - c.y) - (c.y - a.y) * (d.x - c.x)) / det;
    double u = ((c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x)) / det;
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        intersect->x = a.x + t * (b.x - a.x);
        intersect->y = a.y + t * (b.y - a.y);
        return true;
    }
    return false;
}

bool points_are_equal(Point a, Point b) {
    return (fabs(a.x - b.x) < EPSILON) && (fabs(a.y - b.y) < EPSILON);
}