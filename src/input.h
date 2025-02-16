#ifndef INPUT_H
#define INPUT_H

#include "input.h"
#include "point.h"
#include "nearwall.h"

typedef struct {
    char *curve;
    double cell_size;
    int rows;
    int cols;
    Point center;
    double rotation_angle;
    Point rotation_center;
    bool enable_nwl;
    NearWallLayer nwl;
    char *outputfile;
} Input;

Input get_input();

#endif