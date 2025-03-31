#ifndef INPUT_H
#define INPUT_H

#include "input.h"
#include "point.h"
#include "nearwall.h"
#include "shock.h"

typedef struct {
    char *curve;
    char *external_curve;
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
    bool enable_external_nwl;
    NearWallLayer nwl;
    NearWallLayer external_nwl;
    NearShockLayer nsl;
    char *outputfile;
    bool enable_boundaries;
} Input;

void get_input(int argc, char *argv[]);

extern Input input;

#endif