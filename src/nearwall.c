#include "nearwall.h"
#include "point.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846
#define EPSILON 1e-8

double get_SF(NearWallLayer nwl) {
    if (nwl.last > nwl.first) {
        if (nwl.n < 3) return 1.0;
        if (nwl.n == 3) return (nwl.distance-nwl.last-nwl.first)/nwl.first;
        if (nwl.n == 4) return nwl.last/(nwl.distance-nwl.first);
        // TODO: SF(n, L)
        return 1.0;
    }
    if (nwl.n < 3) return 1.0;
    // TODO: SF(n)
    return 1.0;
}

bool compute_offset(Point **offset, int *n_offset, Point *body, int n_body, double nwl_distance) {

    double A = 0.0;
    for (int i = 0; i < n_body; i++) {
        int j = (i + 1) % n_body;
        A += (body[i].x * body[j].y) - (body[j].x * body[i].y);
    }
    int v = 2 * (A < 0) - 1;

    int k = n_body;

    for (int i = 0; i < n_body; i++) {
        int h = (n_body + i - 1) % n_body;
        int j = (i + 1) % n_body;

        double tx_h = (body[i].x - body[h].x);
        double ty_h = (body[i].y - body[h].y);
        double lh = sqrt(tx_h * tx_h + ty_h * ty_h);
        tx_h /= lh;
        ty_h /= lh;

        double tx_j = (body[j].x - body[i].x);
        double ty_j = (body[j].y - body[i].y);
        double lj = sqrt(tx_j * tx_j + ty_j * ty_j);
        tx_j /= lj;
        ty_j /= lj;

        double angle = acos(tx_h * tx_j + ty_h * ty_j);
        double cross_product = tx_h * ty_j - ty_h * tx_j;
        if (cross_product * v > 0) {
            angle = 2 * PI - angle;
        }

        if (angle < EPSILON || angle > 2 * PI - EPSILON) return false;
        if (angle > 1.5 * PI) {
            body[i].x = 0.5 * (body[h].x + body[j].x);
            body[i].y = 0.5 * (body[h].y + body[j].y);
        }
        if (angle < 0.75 * PI) k += 2;
        // TODO: Improve offset quality
    }

    *n_offset = k;
    *offset = (Point *)malloc((*n_offset) * sizeof(Point));

    k = 0;
    for (int i = 0; i < n_body; i++) {
        int h = (n_body + i - 1) % n_body;
        int j = (i + 1) % n_body;

        double tx_h = (body[i].x - body[h].x);
        double ty_h = (body[i].y - body[h].y);
        double lh = sqrt(tx_h * tx_h + ty_h * ty_h);
        tx_h /= lh;
        ty_h /= lh;

        double tx_j = (body[j].x - body[i].x);
        double ty_j = (body[j].y - body[i].y);
        double lj = sqrt(tx_j * tx_j + ty_j * ty_j);
        tx_j /= lj;
        ty_j /= lj;

        double bx = tx_h + tx_j;
        double by = ty_h + ty_j;
        double lb = sqrt(bx * bx + by * by);
        bx /= lb;
        by /= lb;

        double nx = -by * v;
        double ny = bx * v;
    
        double angle = acos(tx_h * tx_j + ty_h * ty_j);
        double cross_product = tx_h * ty_j - ty_h * tx_j;
        if (cross_product * v > 0) {
            angle = 2 * PI - angle;
        }

        if (angle < 0.75 * PI) {
            double n1x = -ty_h * v;
            double n1y = tx_h * v;

            (*offset)[k].x = body[i].x + n1x * nwl_distance;
            (*offset)[k].y = body[i].y + n1y * nwl_distance;
            k++;

            (*offset)[k].x = body[i].x + nx * nwl_distance;
            (*offset)[k].y = body[i].y + ny * nwl_distance;
            k++;

            double n2x = -ty_j * v;
            double n2y = tx_j * v;

            (*offset)[k].x = body[i].x + n2x * nwl_distance;
            (*offset)[k].y = body[i].y + n2y * nwl_distance;
            k++;
        } else if (angle > 1.25 * PI) {
            (*offset)[k].x = body[i].x + sqrt(2) * nx * nwl_distance;
            (*offset)[k].y = body[i].y + sqrt(2) * ny * nwl_distance;
            k++;       
        } else {
            (*offset)[k].x = body[i].x + nx * nwl_distance;
            (*offset)[k].y = body[i].y + ny * nwl_distance;
            k++;
        }
    }
    return true;
}

double get_distance(Node a, Node b) {
    double dx = a.position.x - b.position.x;
    double dy = a.position.y - b.position.y;
    return sqrt(dx*dx + dy*dy);
}

int nearest_node(Node *nodes, int *ids, int n, int current, int *flag, double cell_size) {
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

Point nearest_point(Point p, Point *body, int n_points, double nwl_distance) {
    Point nearest;
    double min_distance = 10 * nwl_distance;
    double ABx, ABy, APx, APy, t, Xq, Yq, dx, dy, distance;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;

        ABx = body[j].x - body[i].x;
        ABy = body[j].y - body[i].y;
        APx = p.x - body[i].x;
        APy = p.y - body[i].y;
        t = (APx * ABx + APy * ABy)/(ABx * ABx + ABy * ABy);
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        Xq = body[i].x + t * ABx;
        Yq = body[i].y + t * ABy;
        dx = Xq - p.x;
        dy = Yq - p.y;
        distance = sqrt(dx * dx + dy * dy);

        if (distance < min_distance) {
            min_distance = distance;
            nearest.x = Xq;
            nearest.y = Yq;
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
        int nearest = nearest_node(type2_nodes, ids, *n_offset_nodes, current, flag, cell_size);
        if (nearest == -1) break;

        flag[nearest] = 1;
        (*offset_nodes)[i] = ids[nearest];
        current = nearest;
    }
}

void extrude_near_wall_cells(
    Element **elements, int *n_elements,
    Node **nodes, int *n_nodes,
    Point *body, int n_body,
    int *offset_nodes, int n_offset_nodes,
    NearWallLayer nwl
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
    bool clockwise = (A < 0);

    Point pB[n_offset_nodes];
    for (int i = 0; i < n_offset_nodes; i++) {
        pB[i] = nearest_point(pA[i], body, n_body, nwl.distance);
    }

    // TODO: Optimize surface points distribution

    if (nwl.n == 1) {
        int kk = *n_nodes;

        for (int i = 0; i < n_offset_nodes; i++) {
            (*nodes)[*n_nodes].position = pB[i];
            (*nodes)[*n_nodes].type = 3;
            (*nodes)[*n_nodes].id = *n_nodes + 1;
            (*n_nodes)++;
        }

        for (int i = 0; i < n_offset_nodes; i++) {
            int n1, n2, n3, n4;
            if (clockwise) {
                n1 = (*nodes)[kk + i].id;
                n2 = (*nodes)[kk + (i + 1) % n_offset_nodes].id;
                n3 = offset_nodes[(i + 1) % n_offset_nodes];
                n4 = offset_nodes[i];
            } else {
                n4 = (*nodes)[kk + i].id;
                n3 = (*nodes)[kk + (i + 1) % n_offset_nodes].id;
                n2 = offset_nodes[(i + 1) % n_offset_nodes];
                n1 = offset_nodes[i];
            }
            (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            (*n_elements)++;
        }
        return;
    }

    // TODO: Create Nodes and Elements with geometric progression distribution
    double x[nwl.n - 1];
    for (int i = 0; i < nwl.n - 1; i++) {
        double K = 0;
        for (int k = 0; k <= i; k++) {
            K += pow(nwl.SF, k);
        }
        x[i] = nwl.first * K;
    }
}