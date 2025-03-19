#ifndef INPUT_H
#define INPUT_H

#include "input.h"
#include "point.h"
#include "nearwall.h"
#include "shock.h"

typedef struct {
    char *curve;
    char *shock_curve;
    double cell_size;
    int coarsening_levels[4];
    int coarsening_cells[4];
    bool fast_coarsening;
    bool conformal_coarsening;
    int rows;
    int cols;
    Point center;
    double rotation_angle;
    Point rotation_center;
    bool smoothing;
    int smoothing_iterations;
    bool enable_nwl;
    NearWallLayer nwl;
    NearShockLayer nsl;
    char *outputfile;
} Input;

Input get_input(int argc, char *argv[]);

#endif