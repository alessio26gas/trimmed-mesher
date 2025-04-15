#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "input.h"
#include "point.h"


#define PI 3.14159265358979323846
#define EPSILON 1e-12

int load_points(const char *filename, Point **points, int *n_points) {

    if (!strcmp(filename, "")) {
        *n_points = 0;
        return 0;
    }

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

    *points = (Point *)malloc(count * sizeof(Point));
    if (!*points) {
        perror("An error has occurred");
        fclose(file);
        return -1;
    }

    rewind(file);
    int i = 0;
    while (fscanf(file, "%lf,%lf", &(*points)[i].x, &(*points)[i].y) == 2) {
        i++;
    }

    *n_points = count;
    fclose(file);
    return 0;
}

void rotate_body()
{
    double AoA = input.rotation_angle / 180 * PI;
    double xc = input.rotation_center.x;
    double yc = input.rotation_center.y;
    double xt, yt;
    for (int i = 0; i < n_body; i++) {
        xt = cos(AoA) * (body[i].x - xc) + sin(AoA) * (body[i].y - yc);
        yt = -sin(AoA) * (body[i].x - xc) + cos(AoA) * (body[i].y - yc);
        body[i].x = xt + xc;
        body[i].y = yt + yc;
    }
}

int get_point_type(Point p, Point *body, int n_points) {
    int inside = 0;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;

        if (points_are_equal(p, body[i])) return 2;

        if (((body[i].y > p.y) != (body[j].y > p.y)) &&
            (p.x < (body[j].x - body[i].x) * (p.y - body[i].y) / (body[j].y - body[i].y) + body[i].x)) {
            inside = !inside;
        }
    }
    return inside;
}

bool is_near_body(Point *p, Point *body, int n_points, double cell_size) {
    double min_distance = cell_size;
    double ABx, ABy, APx, APy, t, Xq, Yq, dx, dy, distance;
    double xp = (*p).x, yp = (*p).y;
    int i_min = 0;
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
            i_min = i;
        }
    }

    double xp2 = xp, yp2 = yp;
    for (int i = 0; i < n_points; i++) {
        if (i == i_min) continue;
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

        if (distance < min_distance + EPSILON) {
            xp2 = Xq;
            yp2 = Yq;
        }
    }

    if (min_distance < cell_size / (2.0 + EPSILON)) {
        (*p).x = (xp + xp2) / 2;
        (*p).y = (yp + yp2) / 2;
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

bool points_are_close(Point a, Point b, double toll) {
    return (fabs(a.x - b.x) < toll + EPSILON) && (fabs(a.y - b.y) < toll + EPSILON);
}

double get_distance(Point a, Point b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx*dx + dy*dy);
}