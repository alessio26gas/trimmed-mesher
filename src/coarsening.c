#include "coarsening.h"
#include <math.h>
#include <stdbool.h>

#define EPSILON 1e-12

void coarsening(Node **nodes, int *n_nodes, Element **elements, int *n_elements, Input input) {

    bool conformal = false;

    int max_level = 0;
    for (int i = 0; i < 4; i++) {
        if (input.coarsening_levels[i] > max_level) {
            max_level = input.coarsening_levels[i];
        }
    }

    double xL = input.center.x - input.cell_size * input.cols / 2;
    double xH = input.center.x + input.cell_size * input.cols / 2;
    double yL = input.center.y - input.cell_size * input.rows / 2;
    double yH = input.center.y + input.cell_size * input.rows / 2;

    for (int l = 1; l <= max_level; l++) {

        double cell_size = input.cell_size * pow(2, l);

        for (int side = 0; side < 4; side++) {

            if (input.coarsening_levels[side] < l) continue;

            int K = (*n_nodes);

            int n = round(((side % 2) * (yH - yL) + ((side + 1) % 2) * (xH - xL)) / cell_size);

            for (int j = conformal; j <= input.coarsening_cells[side] + conformal; j++) {
                for (int i = 0; i <= n; i++) {
                    Point p;
                    if (side==0) p = (Point){xL + i * cell_size, yL - j * cell_size};
                    if (side==1) p = (Point){xH + j * cell_size, yL + i * cell_size};
                    if (side==2) p = (Point){xH - i * cell_size, yH + j * cell_size};
                    if (side==3) p = (Point){xL - j * cell_size, yH - i * cell_size};
                    (*nodes)[*n_nodes] = (Node){*n_nodes + 1, 0, p};
                    (*n_nodes)++;
                }
            }

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < input.coarsening_cells[side]; j++) {
                    (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {
                        K + 1 + i     +  j      * (n + 1),
                        K + 1 + i     + (j + 1) * (n + 1),
                        K + 1 + i + 1 + (j + 1) * (n + 1),
                        K + 1 + i + 1 +  j      * (n + 1)
                    }};
                    (*n_elements)++;
                }
            }

            if (conformal) {
                // TODO: Identify Boundary Nodes

                // TODO: Create Triangular Elements
            }
        }

        double yLold = yL;
        double xHold = xH;
        double yHold = yH;
        double xLold = xL;
        yL -= (input.coarsening_cells[0] + conformal) * cell_size * (input.coarsening_levels[0] >= l);
        xH += (input.coarsening_cells[1] + conformal) * cell_size * (input.coarsening_levels[1] >= l);
        yH += (input.coarsening_cells[2] + conformal) * cell_size * (input.coarsening_levels[2] >= l);
        xL -= (input.coarsening_cells[3] + conformal) * cell_size * (input.coarsening_levels[3] >= l);

        if (yL < yLold - EPSILON) {
            if (xL < xLold - EPSILON) {
                // nodes

                // elements
            }
            if (xH > xHold + EPSILON) {
                // corner in basso a destra
            }
        }
        if (yH > yHold + EPSILON) {
            if (xL < xLold - EPSILON) {
                // corner in alto a sinistra
            }
            if (xH > xHold + EPSILON) {
                // corner in alto a destra
            }
        }

        for (int i = 0; i < 4; i++)
            input.coarsening_cells[i] /= 2;
    }
}