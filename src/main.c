#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "main.h"
#include "input.h"
#include "point.h"
#include "node.h"
#include "element.h"
#include "export.h"
#include "nearwall.h"
#include "shock.h"
#include "smoothing.h"
#include "coarsening.h"

#define PI 3.14159265358979323846
#define EPSILON 1e-12

int n_nodes;
Node *nodes;

int n_elements;
Element *elements;

int n_boundaries;
Element *boundaries;

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

    get_input(argc, argv);

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

    Point *external;
    int n_external;
    if (strcmp(input.external_curve, "")) printf("Loading external curve...");
    start = end;
    if (load_points(input.external_curve, &external, &n_external) != 0) {
        printf(" Failed.\n");
        return -1;
    }
    end = clock();
    if (n_external > 0) printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    bool enable_external_nwl = input.enable_external_nwl;
    NearWallLayer external_nwl;
    Point *external_offset;
    int n_external_offset;
    if (enable_external_nwl) {

        printf("Computing near external boundary parameters...");
        start = end;
        external_nwl = input.external_nwl;
        compute_nwl(&external_nwl);
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

        printf("Computing external boundary offset...");
        start = end;
        if (!compute_offset(&external_offset, &n_external_offset, external, n_external, external_nwl.distance, cell_size)) {
            printf("Invalid geometry.\n");
            return -1;
        }
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    } else {
        external_offset = external;
        n_external_offset = n_external;
    }

    Point *shock;
    int n_shock;
    if (strcmp(input.shock_curve, "")) printf("Loading shockwave curve...");
    start = end;
    if (load_points(input.shock_curve, &shock, &n_shock) != 0) {
        printf(" Failed.\n");
        return -1;
    }
    end = clock();
    if (n_shock > 0) printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    NearShockLayer nsl;
    Point *shock_offset;
    int n_shock_offset;
    if (n_shock > 0) {

        printf("Computing near shockwave parameters...");
        start = end;
        nsl = input.nsl;
        compute_nsl(&nsl);
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

        printf("Computing shockwave offset...");
        start = end;
        if (!compute_shock_offset(&shock_offset, &n_shock_offset, shock, n_shock, nsl.distance, cell_size)) {
            printf("Invalid shockwave geometry.\n");
            return -1;
        }
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    }

    printf("Allocating memory for nodes...");
    start = end;
    n_nodes = (rows + 1) * (cols + 1);
    nodes = malloc(2 * n_nodes * sizeof(Node)); // 2?
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
    if (strcmp(input.external_curve, ""))
    {
        node_id = 0;
        for (int i = 0; i <= rows; i++) {
            for (int j = 0; j <= cols; j++) {
                node_id++;
                if (nodes[node_id - 1].type != 0) continue;
                Point p = nodes[node_id - 1].position;
                int type = get_point_type(p, external_offset, n_external_offset);
                if (type == 0) nodes[node_id - 1].type = 100;
                if (type == 2) nodes[node_id - 1].type = 200;
            }
        }
    }
    if (strcmp(input.shock_curve, ""))
    {
        for (int i = 0; i <= rows; i++) {
            for (int j = 0; j <= cols; j++) {
                node_id++;
                if (nodes[node_id - 1].type != 0) continue;
                Point p = nodes[node_id - 1].position;
                int type = get_point_type(p, shock_offset, n_shock_offset);
                if (type == 0) nodes[node_id - 1].type = 10;
                if (type == 2) nodes[node_id - 1].type = 20;
            }
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
        if (strcmp(input.shock_curve, ""))
        if (is_near_body(&(nodes[i].position), shock_offset, n_shock_offset, cell_size)) {
            nodes[i].type = 20;
        }
        if (strcmp(input.external_curve, ""))
        if (is_near_body(&(nodes[i].position), external_offset, n_external_offset, cell_size)) {
            nodes[i].type = 200;
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
            if (strcmp(input.shock_curve, ""))
            if (is_near_body(&(nodes[i].position), shock_offset, n_shock_offset, cell_size)) {
                nodes[i].type = 20;
            }
            if (strcmp(input.external_curve, ""))
            if (is_near_body(&(nodes[i].position), external_offset, n_external_offset, cell_size)) {
                nodes[i].type = 200;
            }
        }
    }
    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    printf("Allocating memory for elements...");
    start = end;
    elements = malloc(2 * rows * cols * sizeof(Element)); // 2?
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
    int element_id = 1;
    n_elements = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int n1 = i * (cols + 1) + j + 1;
            int n2 = n1 + 1;
            int n3 = n1 + (cols + 1) + 1;
            int n4 = n1 + (cols + 1);

            int ex = (nodes[n1-1].type==0) + (nodes[n2-1].type==0) + (nodes[n3-1].type==0) + (nodes[n4-1].type==0);
            int in = (nodes[n1-1].type==1) + (nodes[n2-1].type==1) + (nodes[n3-1].type==1) + (nodes[n4-1].type==1);
            int s_ex = (nodes[n1-1].type==10) + (nodes[n2-1].type==10) + (nodes[n3-1].type==10) + (nodes[n4-1].type==10);

            if (ex == 0) continue;
            if (s_ex == 4) continue;

            if (in == 0 && s_ex == 0) {
                elements[n_elements++] = (Element){element_id++, 3, 4, {n1, n2, n3, n4}};
                continue;
            }

            if (s_ex == 0) {
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
            } else {
                if (nodes[n1-1].type==10) (nodes[n1-1].type=1);
                if (nodes[n2-1].type==10) (nodes[n2-1].type=1);
                if (nodes[n3-1].type==10) (nodes[n3-1].type=1);
                if (nodes[n4-1].type==10) (nodes[n4-1].type=1);
                if (nodes[n1-1].type==20) (nodes[n1-1].type=2);
                if (nodes[n2-1].type==20) (nodes[n2-1].type=2);
                if (nodes[n3-1].type==20) (nodes[n3-1].type=2);
                if (nodes[n4-1].type==20) (nodes[n4-1].type=2);
                switch (ex) {
                    case 3:
                        pentagons = true;
                        int n5;
                        int flag = penta_vertices(&n1, &n2, &n3, &n4, &n5, shock_offset, n_shock_offset, &nodes, &n_nodes);
                        elements[n_elements++] = (Element){element_id++, 4, 5, {n1, n2, n3, n4, n5}, flag};
                        break;
                    case 2:
                        quad_vertices(&n1, &n2, &n3, &n4, shock_offset, n_shock_offset, &nodes, &n_nodes);
                        elements[n_elements++] = (Element){element_id++, 3, 4, {n1, n2, n3, n4}};
                        break;
                    case 1:
                        tria_vertices(&n1, &n2, &n3, &n4, shock_offset, n_shock_offset, &nodes, &n_nodes);
                        elements[n_elements++] = (Element){element_id++, 2, 3, {n1, n2, n3}};
                        break;
                }
                if (nodes[n1-1].type==1) (nodes[n1-1].type=10);
                if (nodes[n2-1].type==1) (nodes[n2-1].type=10);
                if (nodes[n3-1].type==1) (nodes[n3-1].type=10);
                if (nodes[n4-1].type==1) (nodes[n4-1].type=10);
                if (nodes[n1-1].type==2) (nodes[n1-1].type=20);
                if (nodes[n2-1].type==2) (nodes[n2-1].type=20);
                if (nodes[n3-1].type==2) (nodes[n3-1].type=20);
                if (nodes[n4-1].type==2) (nodes[n4-1].type=20);
            }
        }
    }
    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    if (pentagons) {
        printf("Splitting pentagons...");
        start = end;
        split_pentagons();
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    }

    smoothing();

    coarsening();

    if (enable_nwl) {
        printf("Generating near wall cells...");
        start = clock();
        int *offset_nodes;
        int n_offset_nodes = 0;
        int simmetry = get_offset_nodes(&offset_nodes, &n_offset_nodes, cell_size, X0, Y0, rows, cols);
        extrude_near_wall_cells(
            body, n_body,
            offset_nodes, n_offset_nodes,
            nwl, simmetry,
            X0, cols, cell_size
        );
        free(offset_nodes);
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    }

    if (n_shock > 0) {
        printf("Generating near shock cells...");
        start = end;
        int *offset_nodes;
        int n_offset_nodes = 0;
        for (int i = 0; i < n_nodes; i++) {
            if (nodes[i].type == 2) nodes[i].type = 200;
            if (nodes[i].type == 20) nodes[i].type = 2;
        }
        int simmetry = get_offset_nodes(&offset_nodes, &n_offset_nodes, cell_size, X0, Y0, rows, cols);
        extrude_near_shock_cells(
            &elements, &n_elements,
            &nodes, &n_nodes,
            shock, n_shock,
            offset_nodes, n_offset_nodes,
            nsl, simmetry,
            X0, cols, cell_size
        );
        free(offset_nodes);
        for (int i = 0; i < n_nodes; i++) {
            if (nodes[i].type == 2) nodes[i].type = 20;
            if (nodes[i].type == 200) nodes[i].type = 2;
        }
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    }

    bool enable_boundaries = input.enable_boundaries;
    n_boundaries = 0;

    if (enable_boundaries) {
        printf("Computing boundary elements...");
        start = end;
        boundaries = malloc(4 * (rows + cols) * sizeof(Element)); // 4?
        if (!boundaries) {
            perror("An error has occurred");
            printf(" Failed.\n");
            return -1;
        }
    
        int boundary_id = n_elements + 1;
        for (int c = 0; c < cols; c++) {
            if (nodes[c].type == 0 || nodes[c+1].type == 0) {
                int n1 = nodes[c].id;
                int n2 = nodes[c+1].id;
                boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 1};
                n_boundaries++;
                boundary_id++;
            }
        }
        for (int r = 0; r < rows; r++) {
            if (nodes[r * (cols + 1) + cols].type == 0 || nodes[(r + 1) * (cols + 1) + cols].type == 0) {
                int n1 = nodes[r * (cols + 1) + cols].id;
                int n2 = nodes[(r + 1) * (cols + 1) + cols].id;
                boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 2};
                n_boundaries++;
                boundary_id++;
            }
        }
        for (int c = cols; c > 0; c--) {
            if (nodes[c + rows*(cols + 1)].type == 0 || nodes[c-1 + rows*(cols + 1)].type == 0) {
                int n1 = nodes[c + rows*(cols + 1)].id;
                int n2 = nodes[c-1 + rows*(cols + 1)].id;
                boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 3};
                n_boundaries++;
                boundary_id++;
            }
        }
        for (int r = rows; r > 0; r--) {
            if (nodes[r * (cols + 1)].type == 0 || nodes[(r - 1) * (cols + 1)].type == 0) {
                int n1 = nodes[r * (cols + 1)].id;
                int n2 = nodes[(r - 1) * (cols + 1)].id;
                boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 4};
                n_boundaries++;
                boundary_id++;
            }
        }
    
        int on_body = enable_nwl ? 4 : 2;
        int k = 0;
        for (int i = 0; i < n_nodes; i++) {
            if (nodes[i].type == on_body) k++;
        }
        int *ids_on_body = malloc(k * sizeof(int));
        if (!ids_on_body) {
            perror("An error has occurred");
            return 1;
        }
        k = 0;
        for (int i = 0; i < n_nodes; i++) {
            if (nodes[i].type == on_body) {
                ids_on_body[k] = nodes[i].id;
                k++;
            }
        }
    
        int simm = 0;
        for (int i = 0; i < k; i++) {
            if (
                nodes[ids_on_body[i]-1].position.x < X0 + EPSILON ||
                nodes[ids_on_body[i]-1].position.x > X0 + cols * cell_size - EPSILON ||
                nodes[ids_on_body[i]-1].position.y < Y0 + EPSILON ||
                nodes[ids_on_body[i]-1].position.y > Y0 + rows * cell_size - EPSILON
            ) simm = nodes[ids_on_body[i]-1].id;
        }
        if (simm != 0) {
            int index = 0;
            int tmp = ids_on_body[0];
            for (int i = 0; i < k; i++) {
                if (ids_on_body[i] == simm) {
                    index = i;
                    break;
                }
            }
            ids_on_body[0] = simm;
            ids_on_body[index] = tmp;
        }
    
        int *flag = malloc(k * sizeof(int));
        for (int i = 0; i < k; i++) flag[i] = 0;
        int current = 0;
        flag[current] = 1;
    
        int *ids_sorted = malloc(k * sizeof(int));
        for (int i = 0; i < k; i++) ids_sorted[i] = 0;
        Node *nodes_on_body = malloc(k * sizeof(Node));
        for (int i = 0; i < k; i++) {
            nodes_on_body[i] = nodes[ids_on_body[i]-1];
        }
    
        ids_sorted[0] = ids_on_body[current];
    
        for (int i = 1; i < k; i++) {
            int nearest = nearest_node(elements, n_elements, nodes_on_body, ids_on_body, k, current, flag, cell_size);
            if (nearest == -1) break;
    
            flag[nearest] = 1;
            ids_sorted[i] = ids_on_body[nearest];
            current = nearest;
        }
    
        for (int i = 0; i < k - 1; i++) {
            int n1 = ids_sorted[i];
            int n2 = ids_sorted[i+1];
            boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 5};
            n_boundaries++;
            boundary_id++;
        }
    
        if (simm == 0) {
            int n1 = ids_sorted[k-1];
            int n2 = ids_sorted[0];
            boundaries[n_boundaries] = (Element){boundary_id, 1, 2, {n1, n2}, 5};
            n_boundaries++;
            boundary_id++;
        }
    
        end = clock();
        printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);
    }
    
    printf("Writing output mesh file...");
    start = clock();
    write_mesh_file();
    end = clock();
    printf(" Done. (%.2f seconds)\n", (float) (end - start) / CLOCKS_PER_SEC);

    free(nodes);
    free(elements);
    if (enable_boundaries) free(boundaries);
    if (n_body > 0) free(body);
    if (enable_nwl) free(offset);
    if (n_shock > 0) {
        free(shock);
        free(shock_offset);
    }
    if (n_external > 0) free(external);
    if (enable_external_nwl > 0) free(external_offset);

    for (int i = 0; i < 70; i++) printf("-");
    printf("\nMesh completed in %.2f seconds. Nodes: %d. Cells: %d.\n", (float) (clock() - global_start) / CLOCKS_PER_SEC, n_nodes, n_elements);
    for (int i = 0; i < 70; i++) printf("-");
    printf("\n");

    return 0;
}