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
