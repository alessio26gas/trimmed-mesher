#ifndef SHOCK_H
#define SHOCK_H

#include "element.h"

typedef struct {
    double distance;
    double first;
    double last;
    double SF;
    int n;
    int distribution;
} NearShockLayer;

#endif