#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "main.h"
#include "input.h"
#include "point.h"
#include "nearwall.h"
#include "shock.h"

#define PI 3.14159265358979323846
#define EPSILON 1e-12

void compute_nsl(NearShockLayer *nsl) {

    NearWallLayer nwl = {
        .first = nsl->first,
        .last = nsl->last,
        .distance = nsl->distance,
        .n = nsl->n,
        .SF = nsl->SF,
        .min_surf_distance = 0.0,
        .surf_max_iter = 1000,
        .distribution = nsl->distribution,
    };

    compute_nwl(&nwl);
    nsl->first = nwl.first;
    nsl->last = nwl.last;
    nsl->distance = nwl.distance;
    nsl->n = nwl.n;
    nsl->SF = nwl.SF;
}

bool compute_shock_offset(Point **offset, int *n_offset, Point *shock, int n_shock, double nsl_distance, double cell_size) {

    double A = 0.0;
    for (int i = 0; i < n_shock; i++) {
        int j = (i + 1) % n_shock;
        A += (shock[i].x * shock[j].y) - (shock[j].x * shock[i].y);
    }
    int v = 2 * (A > 0) - 1;

    int k = n_shock;

    for (int i = 0; i < n_shock; i++) {
        int h = (n_shock + i - 1) % n_shock;
        int j = (i + 1) % n_shock;

        double tx_h = (shock[h].x - shock[i].x);
        double ty_h = (shock[h].y - shock[i].y);
        double lh = sqrt(tx_h * tx_h + ty_h * ty_h);
        tx_h /= lh;
        ty_h /= lh;

        double tx_j = (shock[j].x - shock[i].x);
        double ty_j = (shock[j].y - shock[i].y);
        double lj = sqrt(tx_j * tx_j + ty_j * ty_j);
        tx_j /= lj;
        ty_j /= lj;

        double angle = acos(tx_h * tx_j + ty_h * ty_j);
        double cross_product = tx_h * ty_j - ty_h * tx_j;
        angle += (cross_product * v < 0) * (2*PI - 2*angle);

        if (angle < EPSILON || angle > 2 * PI - EPSILON) return false;
        if (angle < PI) k += 2 * ceil((PI - angle) * nsl_distance / cell_size);
    }

    *n_offset = k;
    *offset = (Point *)malloc((*n_offset) * sizeof(Point));

    k = 0;
    for (int i = 0; i < n_shock; i++) {
        int h = (n_shock + i - 1) % n_shock;
        int j = (i + 1) % n_shock;

        double tx_h = (shock[h].x - shock[i].x);
        double ty_h = (shock[h].y - shock[i].y);
        double lh = sqrt(tx_h * tx_h + ty_h * ty_h);
        tx_h /= lh;
        ty_h /= lh;

        double tx_j = (shock[j].x - shock[i].x);
        double ty_j = (shock[j].y - shock[i].y);
        double lj = sqrt(tx_j * tx_j + ty_j * ty_j);
        tx_j /= lj;
        ty_j /= lj;

        double angle = acos(tx_h * tx_j + ty_h * ty_j);
        double cross_product = tx_h * ty_j - ty_h * tx_j;
        angle += (cross_product * v < 0) * (2*PI - 2*angle);

        if (angle < PI) {
            double nhx = ty_h * v;
            double nhy = -tx_h * v;

            double njx = -ty_j * v;
            double njy = tx_j * v;

            double phih = atan2(nhy, nhx);
            double phij = atan2(njy, njx);

            int N = 1 + 2 * ceil((PI - angle) * nsl_distance / cell_size);

            for (int n = 0; n < N; n++) {
                double alpha = ((double) n) / (N - 1);
                double phin = atan2(
                    (1 - alpha) * sin(phih) + alpha * sin(phij),
                    (1 - alpha) * cos(phih) + alpha * cos(phij)
                );
                (*offset)[k].x = shock[i].x + cos(phin) * nsl_distance;
                (*offset)[k].y = shock[i].y + sin(phin) * nsl_distance;
                k++;
            }

        } else {
            double bx = tx_h + tx_j;
            double by = ty_h + ty_j;
            double lb = sqrt(bx * bx + by * by);
    
            if (lb > EPSILON) {
                bx /= lb;
                by /= lb;    
            } else {
                bx = -ty_j * v;
                by = tx_j * v;
            }

            (*offset)[k].x = shock[i].x + bx * nsl_distance / sin(angle/2);
            (*offset)[k].y = shock[i].y + by * nsl_distance / sin(angle/2);
            k++;
        }
    }
    return true;
}

void extrude_near_shock_cells(
    Point *shock, int n_shock,
    int *offset_nodes, int n_offset_nodes,
    NearShockLayer nsl, int simmetry,
    double X0, int cols, double cell_size
) {
    Point *pA = malloc(n_offset_nodes * sizeof(Point));
    for (int i = 0; i < n_offset_nodes; i++) {
        pA[i] = nodes[offset_nodes[i]-1].position;
    }

    Point *pB = malloc(n_offset_nodes * sizeof(Point));
    for (int i = 0; i < n_offset_nodes; i++) {
        pB[i] = nearest_point(pA[i], shock, n_shock, nsl.distance);
    }

    int *r = malloc(n_offset_nodes * sizeof(int));
    int *l = malloc(n_offset_nodes * sizeof(int));
    int *c = malloc(n_offset_nodes * sizeof(int));
    for (int i = 0; i < n_offset_nodes; i++) {
        c[i] = 0;
    }
    Point pBold[n_offset_nodes];
    bool singularities = true;

    int n_el = n_offset_nodes;
    int n_el0 = 0;
    if (simmetry == 1) {
        n_el0 = 1;
        n_el = n_offset_nodes - 1;
    }

    for (int i = 0; i < n_offset_nodes; i++) {
        Point pBtemp;
        do {
            pBtemp = pB[i];
            pB[i] = nearest_point(pB[i], shock, n_shock, nsl.distance);
        } while (!points_are_equal(pB[i], pBtemp));
    }

    double A = 0.0;
    for (int i = 0; i < n_offset_nodes; i++) {
        int j = (i + 1) % n_offset_nodes;
        A += (pA[i].x * pA[j].y) - (pA[j].x * pA[i].y);
    }
    bool clockwise = (A > 0); // REVERSED

    if (nsl.n == 1) {
        int kk = n_nodes;

        for (int i = 0; i < n_offset_nodes; i++) {
            nodes[n_nodes].position = pB[i];
            nodes[n_nodes].type = 3;
            nodes[n_nodes].id = n_nodes + 1;
            n_nodes++;
        }

        for (int i = 0; i < n_el; i++) {
            int n1 = nodes[kk + i].id;
            int n2 = nodes[kk + (i + 1) % n_offset_nodes].id;
            int n3 = offset_nodes[(i + 1) % n_offset_nodes];
            int n4 = offset_nodes[i];

            if (clockwise) {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            n_elements++;
        }
        return;
    }

    if (nsl.n == 2) {
        double x = nsl.first / nsl.distance;
        int kk = n_nodes;

        for (int i = 0; i < n_offset_nodes; i++) {
            nodes[n_nodes].position = pB[i];
            nodes[n_nodes].type = 3;
            nodes[n_nodes].id = n_nodes + 1;
            n_nodes++;
        }

        for (int i = 0; i < n_offset_nodes; i++) {
            nodes[n_nodes].position.x = x*pA[i].x + (1-x)*pB[i].x;
            nodes[n_nodes].position.y = x*pA[i].y + (1-x)*pB[i].y;
            nodes[n_nodes].type = 3;
            nodes[n_nodes].id = n_nodes + 1;
            n_nodes++;
        }

        for (int i = 0; i < n_el; i++) {
            int n1 = nodes[kk + i].id;
            int n2 = nodes[kk + (i + 1) % n_offset_nodes].id;
            int n3 = nodes[kk + n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n4 = nodes[kk + n_offset_nodes + i].id;

            if (clockwise) {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            n_elements++;
        }

        for (int i = 0; i < n_el; i++) {
            int n1 = nodes[kk + n_offset_nodes + i].id;
            int n2 = nodes[kk + n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n3 = offset_nodes[(i + 1) % n_offset_nodes];
            int n4 = offset_nodes[i];

            if (clockwise) {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            n_elements++;
        }
        return;
    }

    if (nsl.n == 3 && nsl.last > 0) {
        double x[2];
        x[0] = nsl.first / nsl.distance;
        x[1] = 1 - nsl.last / nsl.distance;
        int kk = n_nodes;

        for (int i = 0; i < n_offset_nodes; i++) {
            nodes[n_nodes].position = pB[i];
            nodes[n_nodes].type = 3;
            nodes[n_nodes].id = n_nodes + 1;
            n_nodes++;
        }

        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < n_offset_nodes; i++) {
                nodes[n_nodes].position.x = x[j]*pA[i].x + (1-x[j])*pB[i].x;
                nodes[n_nodes].position.y = x[j]*pA[i].y + (1-x[j])*pB[i].y;
                nodes[n_nodes].type = 3;
                nodes[n_nodes].id = n_nodes + 1;
                n_nodes++;
            }
        }

        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < n_el; i++) {
                int n1 = nodes[kk + j * n_offset_nodes + i].id;
                int n2 = nodes[kk + j * n_offset_nodes + (i + 1) % n_offset_nodes].id;
                int n3 = nodes[kk + (j + 1) * n_offset_nodes + (i + 1) % n_offset_nodes].id;
                int n4 = nodes[kk + (j + 1) * n_offset_nodes + i].id;

                if (clockwise) {
                    elements[n_elements] = (Element){n_elements + 1, 3, 4, {n1, n2, n3, n4}};
                } else {
                    elements[n_elements] = (Element){n_elements + 1, 3, 4, {n4, n3, n2, n1}};
                }
                n_elements++;
            }
        }

        for (int i = 0; i < n_el; i++) {
            int n1 = nodes[kk + 2 * n_offset_nodes + i].id;
            int n2 = nodes[kk + 2 * n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n3 = offset_nodes[(i + 1) % n_offset_nodes];
            int n4 = offset_nodes[i];

            if (clockwise) {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            n_elements++;
        }
        return;
    }

    double *x = malloc((nsl.n - 1) * sizeof(double));

    if (nsl.distribution == 0) {
        double K = 0;
        for (int i = 0; i < nsl.n - 1; i++) {
            K += pow(nsl.SF, i);
            x[i] = nsl.first / nsl.distance * K;
        }
    } else {
        for (int i = 0; i < nsl.n - 1; i++) {
            x[i] = 1 + tanh(nsl.SF * ((i + 1.0)/(nsl.n) - 1))/tanh(nsl.SF);
        }
    }

    int kk = n_nodes;

    for (int i = 0; i < n_offset_nodes; i++) {
        nodes[n_nodes].position.x = (2*pB[i].x - pA[i].x);
        nodes[n_nodes].position.y = (2*pB[i].y - pA[i].y);
        nodes[n_nodes].type = 3;
        nodes[n_nodes].id = n_nodes + 1;
        n_nodes++;
    }

    for (int j = 0; j < nsl.n - 1; j++) {
        for (int i = 0; i < n_offset_nodes; i++) {
            nodes[n_nodes].position.x = x[j]*pA[i].x + (1-x[j])*(2*pB[i].x - pA[i].x);
            nodes[n_nodes].position.y = x[j]*pA[i].y + (1-x[j])*(2*pB[i].y - pA[i].y);
            nodes[n_nodes].type = 3;
            nodes[n_nodes].id = n_nodes + 1;
            n_nodes++;
        }
    }

    for (int j = 0; j < nsl.n - 1; j++) {
        for (int i = 0; i < n_el; i++) {
            int n1 = nodes[kk + j * n_offset_nodes + i].id;
            int n2 = nodes[kk + j * n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n3 = nodes[kk + (j + 1) * n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n4 = nodes[kk + (j + 1) * n_offset_nodes + i].id;

            if (clockwise) {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                elements[n_elements] = (Element){n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            n_elements++;
        }
    }

    for (int i = 0; i < n_el; i++) {
        int n1 = nodes[kk + (nsl.n - 1) * n_offset_nodes + i].id;
        int n2 = nodes[kk + (nsl.n - 1) * n_offset_nodes + (i + 1) % n_offset_nodes].id;
        int n3 = offset_nodes[(i + 1) % n_offset_nodes];
        int n4 = offset_nodes[i];

        if (clockwise) {
            elements[n_elements] = (Element){n_elements + 1, 3, 4, {n1, n2, n3, n4}};
        } else {
            elements[n_elements] = (Element){n_elements + 1, 3, 4, {n4, n3, n2, n1}};
        }
        n_elements++;
    }
}