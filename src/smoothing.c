#include "main.h"
#include "input.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define EPSILON 1e-12

void smoothing()
{
    double X0 = - input.cell_size * input.cols / 2 + input.center.x;
    double Y0 = - input.cell_size * input.rows / 2 + input.center.y;

    if (input.smoothing && input.smoothing_iterations > 0) {
        printf("Executing smoothing algorithm...\n");
        clock_t start = clock();
        clock_t end;
        for (int k = 1; k <= input.smoothing_iterations; k++) {
            double sum = 0.0;
            for (int c = 1; c < input.cols; c++) {
                for (int r = 1; r < input.rows; r++) {
                    int i = c + r * (input.cols + 1);
                    if (nodes[i].type != 0) continue;
                    sum += fabs((nodes[i+1].position.x + nodes[i-1].position.x + nodes[i-input.cols-1].position.x + nodes[i+input.cols+1].position.x)/4.0 - nodes[i].position.x);
                    sum += fabs((nodes[i+1].position.y + nodes[i-1].position.y + nodes[i-input.cols-1].position.y + nodes[i+input.cols+1].position.y)/4.0 - nodes[i].position.y);
                    nodes[i].position.x = (nodes[i+1].position.x + nodes[i-1].position.x + nodes[i-input.cols-1].position.x + nodes[i+input.cols+1].position.x)/4.0;
                    nodes[i].position.y = (nodes[i+1].position.y + nodes[i-1].position.y + nodes[i-input.cols-1].position.y + nodes[i+input.cols+1].position.y)/4.0;    
                }
            }
            for (int c = 1; c < input.cols; c++) {
                if (nodes[c].type == 0 && input.coarsening_levels[0] == 0) {
                    nodes[c].position.x = (nodes[c+1].position.x + nodes[c-1].position.x)/2.0;
                }
                int i = c + input.rows*(input.cols + 1);
                if (nodes[i].type == 0 && input.coarsening_levels[2] == 0) {
                    nodes[i].position.x = (nodes[i+1].position.x + nodes[i-1].position.x)/2.0;
                } 
            }
            for (int r = 1; r < input.rows; r++) {
                int i = r * (input.cols + 1);
                if (nodes[i].type == 0 && input.coarsening_levels[1] == 0) {
                    nodes[i].position.y = (nodes[i-input.cols-1].position.y + nodes[i+input.cols+1].position.y)/2.0;
                }
                i = r * (input.cols + 1) + input.cols;
                if (nodes[i].type == 0 && input.coarsening_levels[3] == 0) {
                    nodes[i].position.y = (nodes[i-input.cols-1].position.y + nodes[i+input.cols+1].position.y)/2.0;
                }
            }
            if (sum < EPSILON) {
                printf("Smoothing iteration %d: residuals = %.4e\n", k, sum);
                printf("Smoothing algorithm converged in %d iterations. (%.2f seconds)\n", k, (float) (clock() - start) / CLOCKS_PER_SEC);
                break;
            }
            if (k % 1000 == 0) printf("Smoothing iteration %d: residuals = %.4e\n", k, sum);
            if (k == input.smoothing_iterations) {
                printf("Smoothing algorithm stopped at %d iterations. (%.2f seconds)\n", k, (float) (clock() - start) / CLOCKS_PER_SEC);
            }
        }
    }

    for (int c = 0; c < input.cols + 1; c++) {
        nodes[c].position.y = Y0;
        nodes[c + input.rows*(input.cols + 1)].position.y = Y0 + input.rows * input.cell_size;
    }

    for (int r = 0; r < input.rows + 1; r++) {
        nodes[r * (input.cols + 1)].position.x = X0;
        nodes[r * (input.cols + 1) + input.cols].position.x = X0 + input.cols * input.cell_size;
    }
}