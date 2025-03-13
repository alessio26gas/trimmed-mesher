#ifndef SHOCK_H
#define SHOCK_H

#include "nearwall.h"

typedef struct {
    double distance;
    double first;
    double last;
    double SF;
    int n;
    int distribution;
} NearShockLayer;

void compute_nsl(NearShockLayer *nsl);

#endif