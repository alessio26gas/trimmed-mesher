#include "input.h"
#include <stdio.h>
#include <stdlib.h>

Input get_input() {
    Input input = (Input){

        .curve = "curve.csv",

        .cell_size = 0.01,

        .coarsening_levels = {0, 0, 0, 0},
        .coarsening_cells = {0, 0, 0, 0},
        .fast_coarsening = false,
        .conformal_coarsening = false,

        .rows = 256,
        .cols = 256,

        .center = (Point){0, 0},

        .rotation_angle = 0.0,
        .rotation_center = (Point){0.0, 0.0},

        .enable_nwl = true,
        .nwl.first = 5.0e-5,
        .nwl.last = input.cell_size,
        .nwl.distance = 0,
        .nwl.n = 0,
        .nwl.SF = 1.3,

        .outputfile = "mesh.msh"
    };

    return input;
}