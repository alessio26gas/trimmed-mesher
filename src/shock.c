#include "shock.h"
#include "nearwall.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
