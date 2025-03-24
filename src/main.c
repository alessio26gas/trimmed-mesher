#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "point.h"
#include "node.h"
#include "element.h"
#include "export.h"
#include "nearwall.h"
#include "coarsening.h"
#include "input.h"

#define PI 3.14159265358979323846
#define EPSILON 1e-12

int main(int argc, char *argv[]) {

    setbuf(stdout, NULL);

    for (int i = 0; i < 70; i++) printf("-");
    printf("\n");
    for (int i = 0; i < 17; i++) printf(" ");
    printf("[Trimmed Mesher Development Version]\n");
    for (int i = 0; i < 70; i++) printf("-");
    printf("\n");

    clock_t global_start = clock();
    clock_t start = global_start;

    Input input = get_input(argc, argv);

    double cell_size = input.cell_size;
    int rows = input.rows;
    int cols = input.cols;

    if (rows < EPSILON || cols < EPSILON || cell_size < EPSILON) {
        printf("Invalid input parameters.\n");
        return -1;
    }

    double X0 = - cell_size * cols / 2 + input.center.x;
    double Y0 = - cell_size * rows / 2 + input.center.y;

    Point *body;
    int n_body;
    if (strcmp(input.curve, "")) printf("Loading curve...");
    if (load_points(input.curve, &body, &n_body) != 0) {
        printf(" Failed.\n");
        return -1;
    }
    clock_t end = clock();
    if (n_body > 0) printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    double AoA = input.rotation_angle / 180 * PI;
    double xc = input.rotation_center.x;
    double yc = input.rotation_center.y;
    if (AoA != 0.0) {
        printf("Rotating curve...");
        start = end;
        double xt, yt;
        for (int i = 0; i < n_body; i++) {
            xt = cos(AoA) * (body[i].x - xc) + sin(AoA) * (body[i].y - yc);
            yt = -sin(AoA) * (body[i].x - xc) + cos(AoA) * (body[i].y - yc);
            body[i].x = xt + xc;
            body[i].y = yt + yc;
        }
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    }

    bool enable_nwl = input.enable_nwl;
    NearWallLayer nwl;
    Point *offset;
    int n_offset;
    if (enable_nwl) {

        printf("Computing near wall parameters...");
        start = end;
        nwl = input.nwl;
        compute_nwl(&nwl);
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

        printf("Computing curve offset...");
        start = end;
        if (!compute_offset(&offset, &n_offset, body, n_body, nwl.distance, cell_size)) {
            printf("Invalid geometry.\n");
            return -1;
        }
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    } else {
        offset = body;
        n_offset = n_body;
    }

    printf("Allocating memory for nodes...");
    start = end;
    int n_nodes = (rows + 1) * (cols + 1);
    Node *nodes = malloc(2 * n_nodes * sizeof(Node)); // 2?
    if (!nodes) {
        perror("An error has occurred");
        printf(" Failed.\n");
        return -1;
    }
    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    printf("Generating initial cartesian grid...");
    start = end;
    int node_id = 1;
    for (int i = 0; i <= rows; i++) {
        for (int j = 0; j <= cols; j++) {
            Point p = {X0 + j * cell_size, Y0 + i * cell_size};
            int type = get_point_type(p, offset, n_offset);
            nodes[node_id - 1] = (Node){node_id, type, p};
            node_id++;
        }
    }
    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    printf("Optimizing near body nodes position...");
    start = end;
    for (int i = 0; i < n_nodes; i++) {
        if (nodes[i].type != 0) continue;
        if (is_near_body(&(nodes[i].position), offset, n_offset, cell_size)) {
            nodes[i].type = 2;
        }
    }

    for (int i = cols+2; i < n_nodes-cols-2; i++) {
        if (nodes[i].type != 1) continue;

        int j1 = nodes[i-1].type + nodes[i-cols-1].type + nodes[i-cols-2].type;
        int j2 = nodes[i+1].type + nodes[i-cols-1].type + nodes[i-cols].type;
        int j3 = nodes[i+1].type + nodes[i+cols+1].type + nodes[i+cols+2].type;
        int j4 = nodes[i-1].type + nodes[i+cols+1].type + nodes[i+cols].type;

        if (j1 == 0 || j2 == 0 || j3 == 0 || j4 == 0) {
            if (is_near_body(&(nodes[i].position), offset, n_offset, cell_size)) {
                nodes[i].type = 2;
            }
        }
    }
    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    printf("Allocating memory for elements...");
    start = end;
    Element *elements = malloc(2 * rows * cols * sizeof(Element)); // 2?
    if (!elements) {
        perror("An error has occurred");
        printf(" Failed.\n");
        return -1;
    }
    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    printf("Computing elements...");
    start = end;
    bool pentagons = false;
    int element_id = 1, n_elements = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int n1 = i * (cols + 1) + j + 1;
            int n2 = n1 + 1;
            int n3 = n1 + (cols + 1) + 1;
            int n4 = n1 + (cols + 1);

            int ex = (nodes[n1-1].type==0) + (nodes[n2-1].type==0) + (nodes[n3-1].type==0) + (nodes[n4-1].type==0);
            int in = (nodes[n1-1].type==1) + (nodes[n2-1].type==1) + (nodes[n3-1].type==1) + (nodes[n4-1].type==1);

            if (ex == 0) continue;

            if (in == 0) {
                elements[n_elements++] = (Element){element_id++, 3, 4, {n1, n2, n3, n4}};
                continue;
            }

            switch (ex) {
                case 3:
                    pentagons = true;
                    int n5;
                    int flag = penta_vertices(&n1, &n2, &n3, &n4, &n5, offset, n_offset, &nodes, &n_nodes);
                    elements[n_elements++] = (Element){element_id++, 4, 5, {n1, n2, n3, n4, n5}, flag};
                    break;
                case 2:
                    quad_vertices(&n1, &n2, &n3, &n4, offset, n_offset, &nodes, &n_nodes);
                    elements[n_elements++] = (Element){element_id++, 3, 4, {n1, n2, n3, n4}};
                    break;
                case 1:
                    tria_vertices(&n1, &n2, &n3, &n4, offset, n_offset, &nodes, &n_nodes);
                    elements[n_elements++] = (Element){element_id++, 2, 3, {n1, n2, n3}};
                    break;
            }
        }
    }
    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    if (pentagons) {
        printf("Splitting pentagons...");
        start = end;
        split_pentagons(&elements, &n_elements);
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    }

    if (input.smoothing && input.smoothing_iterations > 0) {
        printf("Executing smoothing algorithm...\n");
        start = end;
        for (int k = 1; k <= input.smoothing_iterations; k++) {
            double sum = 0.0;
            for (int c = 1; c < cols; c++) {
                for (int r = 1; r < rows; r++) {
                    int i = c + r * (cols + 1);
                    if (nodes[i].type != 0) continue;
                    sum += fabs((nodes[i+1].position.x + nodes[i-1].position.x + nodes[i-cols-1].position.x + nodes[i+cols+1].position.x)/4.0 - nodes[i].position.x);
                    sum += fabs((nodes[i+1].position.y + nodes[i-1].position.y + nodes[i-cols-1].position.y + nodes[i+cols+1].position.y)/4.0 - nodes[i].position.y);
                    nodes[i].position.x = (nodes[i+1].position.x + nodes[i-1].position.x + nodes[i-cols-1].position.x + nodes[i+cols+1].position.x)/4.0;
                    nodes[i].position.y = (nodes[i+1].position.y + nodes[i-1].position.y + nodes[i-cols-1].position.y + nodes[i+cols+1].position.y)/4.0;    
                }
            }
            for (int c = 1; c < cols; c++) {
                if (nodes[c].type == 0 && input.coarsening_levels[0] == 0) {
                    nodes[c].position.x = (nodes[c+1].position.x + nodes[c-1].position.x)/2.0;
                }
                int i = c + rows*(cols + 1);
                if (nodes[i].type == 0 && input.coarsening_levels[2] == 0) {
                    nodes[i].position.x = (nodes[i+1].position.x + nodes[i-1].position.x)/2.0;
                } 
            }
            for (int r = 1; r < rows; r++) {
                int i = r * (cols + 1);
                if (nodes[i].type == 0 && input.coarsening_levels[1] == 0) {
                    nodes[i].position.y = (nodes[i-cols-1].position.y + nodes[i+cols+1].position.y)/2.0;
                }
                i = r * (cols + 1) + cols;
                if (nodes[i].type == 0 && input.coarsening_levels[3] == 0) {
                    nodes[i].position.y = (nodes[i-cols-1].position.y + nodes[i+cols+1].position.y)/2.0;
                }
            }
            if (sum < EPSILON) {
                printf("Smoothing iteration %d: residuals = %.4e\n", k, sum);
                end = clock();
                printf("Smoothing algorithm converged in %d iterations. (%.2f seconds)\n", k, (float) (end - start) / CLOCKS_PER_SEC);
                break;
            }
            if (k % 1000 == 0) printf("Smoothing iteration %d: residuals = %.4e\n", k, sum);
            if (k == input.smoothing_iterations) {
                end = clock();
                printf("Smoothing algorithm stopped at %d iterations. (%.2f seconds)\n", k, (float) (end - start) / CLOCKS_PER_SEC);
            }
        }
    }

    for (int c = 0; c < cols + 1; c++) {
        nodes[c].position.y = Y0;
        nodes[c + rows*(cols + 1)].position.y = Y0 + rows * cell_size;
    }

    for (int r = 0; r < rows + 1; r++) {
        nodes[r * (cols + 1)].position.x = X0;
        nodes[r * (cols + 1) + cols].position.x = X0 + cols * cell_size;
    }

    coarsening(&nodes, &n_nodes, &elements, &n_elements, input);

    if (enable_nwl) {
        printf("Generating near wall cells...");
        start = end;
        int *offset_nodes;
        int n_offset_nodes = 0;
        int simmetry = get_offset_nodes(&offset_nodes, &n_offset_nodes, nodes, n_nodes, elements, n_elements, cell_size, X0, Y0, rows, cols);
        extrude_near_wall_cells(
            &elements, &n_elements,
            &nodes, &n_nodes,
            body, n_body,
            offset_nodes, n_offset_nodes,
            nwl, simmetry,
            X0, cols, cell_size
        );
        free(offset_nodes);
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    }

    printf("Computing boundary elements...");
    start = end;

    Element *boundaries = malloc(4 * (rows + cols) * sizeof(Element)); // 4?
    if (!boundaries) {
        perror("An error has occurred");
        printf(" Failed.\n");
        return -1;
    }

    int n_boundaries = 0;
    int boundary_id = n_elements + 1;
    for (int c = 0; c < cols; c++) {
        if (nodes[c].type != 1 && nodes[c+1].type != 1) {
            int n1 = nodes[c].id;
            int n2 = nodes[c+1].id;
            boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 4};
            n_boundaries++;
            boundary_id++;
        }
    }
    for (int r = 0; r < rows; r++) {
        if (nodes[r * (cols + 1) + cols].type != 1 && nodes[(r + 1) * (cols + 1) + cols].type != 1) {
            int n1 = nodes[r * (cols + 1) + cols].id;
            int n2 = nodes[(r + 1) * (cols + 1) + cols].id;
            boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 5};
            n_boundaries++;
            boundary_id++;
        }
    }
    for (int c = cols; c > 0; c--) {
        if (nodes[c + rows*(cols + 1)].type != 1 && nodes[c-1 + rows*(cols + 1)].type != 1) {
            int n1 = nodes[c + rows*(cols + 1)].id;
            int n2 = nodes[c-1 + rows*(cols + 1)].id;
            boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 6};
            n_boundaries++;
            boundary_id++;
        }
    }
    for (int r = rows; r > 0; r--) {
        if (nodes[r * (cols + 1)].type != 1 && nodes[(r - 1) * (cols + 1)].type != 1) {
            int n1 = nodes[r * (cols + 1)].id;
            int n2 = nodes[(r - 1) * (cols + 1)].id;
            boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 7};
            n_boundaries++;
            boundary_id++;
        }
    }

    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    printf("Writing output mesh file...");
    start = clock();
    write_mesh_file(input.outputfile, nodes, n_nodes, elements, n_elements, boundaries, n_boundaries);
    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    free(nodes);
    free(elements);
    free(boundaries);
    if (n_body > 0) free(body);
    if (enable_nwl) free(offset);

    for (int i = 0; i < 70; i++) printf("-");
    printf("\nMesh completed in %.2f seconds. Nodes: %d. Cells: %d.\n", (float) (clock() - global_start) / CLOCKS_PER_SEC, n_nodes, n_elements);
    for (int i = 0; i < 70; i++) printf("-");
    printf("\n");

    return 0;
}