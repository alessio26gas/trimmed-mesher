#include "input.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Input get_input(int argc, char *argv[]) {

    Input input = {};

    char *curve = "curve.csv";
    char *outputfile = "mesh.msh";

    double cell_size = 0.01;

    int coarsening_levels[] = {0, 0, 0, 0};
    int coarsening_cells[] = {0, 0, 0, 0};
    bool fast_coarsening = false;
    bool conformal_coarsening = false;

    int rows = 256;
    int cols = 256;

    Point center = (Point){0.0, 0};

    double rotation_angle = 0.0;
    Point rotation_center = (Point){0.0, 0.0};

    bool smoothing = true;
    int smoothing_iterations = 3000;

    bool enable_nwl = true;
    NearWallLayer nwl = {
        .first = 5.0e-5,
        .last = cell_size,
        .distance = 0,
        .n = 0,
        .SF = 1.25,
        .min_surf_distance = 0.0,
        .surf_max_iter = 1000,
        .distribution = 0,
    };

    char *shock_curve = "";
    NearShockLayer nsl = {
        .first = 0.00125,
        .last = 0,
        .distance = 0.05,
        .n = 40,
        .SF = 0,
        .distribution = 0,
    };

    if (argc == 32) {
        curve = argv[1];
        outputfile = argv[2];
        cell_size = atof(argv[3]);
        for (int i = 0; i < 4; i++) {
            coarsening_levels[i] = atoi(argv[4 + i]);
            coarsening_cells[i] = atoi(argv[8 + i]);
        }
        fast_coarsening = atoi(argv[12]);
        conformal_coarsening = atoi(argv[13]);
        rows = atoi(argv[14]);
        cols = atoi(argv[15]);
        center.x = atof(argv[16]);
        center.y = atof(argv[17]);
        rotation_angle = atof(argv[18]);
        rotation_center.x = atof(argv[19]);
        rotation_center.y = atof(argv[20]);
        smoothing = atoi(argv[21]);
        smoothing_iterations = atoi(argv[22]);
        enable_nwl = atoi(argv[23]);
        nwl.first = atof(argv[24]);
        nwl.last = atof(argv[25]);
        nwl.distance = atof(argv[26]);
        nwl.n = atoi(argv[27]);
        nwl.SF = atof(argv[28]);
        nwl.min_surf_distance = atof(argv[29]);
        nwl.surf_max_iter = atof(argv[30]);
        nwl.distribution = atoi(argv[31]);
    }

    input = (Input){
        .curve = curve,
        .cell_size = cell_size,
        .fast_coarsening = fast_coarsening,
        .conformal_coarsening = conformal_coarsening,
        .rows = rows,
        .cols = cols,
        .center = center,
        .rotation_angle = rotation_angle,
        .rotation_center = rotation_center,
        .smoothing = smoothing,
        .smoothing_iterations = smoothing_iterations,
        .enable_nwl = enable_nwl,
        .nwl = nwl,
        .outputfile = outputfile,
        .shock_curve = shock_curve,
        .nsl = nsl
    };

    memcpy(input.coarsening_levels, coarsening_levels, sizeof(coarsening_levels));
    memcpy(input.coarsening_cells, coarsening_cells, sizeof(coarsening_cells));

    return input;
}