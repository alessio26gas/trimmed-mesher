#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "point.h"
#include "node.h"
#include "element.h"
#include "export.h"
#include "nearwall.h"

#define EPSILON 1e-8

int main(int argc, char *argv[]) {

    if (argc < 4 || argc > 7) {
        printf("Usage: %s <curve.csv> <cell_size> <rows> [columns] [X0] [Y0]\n", argv[0]);
        return -1;
    }

    double cell_size = atof(argv[2]);
    int rows = atoi(argv[3]);
    int cols;

    if (argc == 4) {
        cols = rows;
    } else {
        cols = atoi(argv[4]);
    }

    if (rows < EPSILON || cols < EPSILON || cell_size < EPSILON) {
        printf("Unvalid values.\n");
        printf("Usage: %s <curve.csv> <cell_size> <rows> [columns] [X0] [Y0]\n", argv[0]);
        return -1;
    }

    double X0 = - cell_size * cols / 2;
    double Y0 = - cell_size * rows / 2;
    if (argc >= 6) X0 += atof(argv[5]);
    if (argc == 7) Y0 += atof(argv[6]);

    Point *body;
    int n_body;
    if (load_body(argv[1], &body, &n_body) != 0) {
        return -1;
    }

    NearWallLayer nwl;
    nwl.first = 0.001;
    nwl.last = 0.8 * cell_size;
    nwl.distance = 1.8 * cell_size;
    nwl.n = 1;
    nwl.SF = get_SF(nwl);

    Point *offset;
    int n_offset;
    if (!compute_offset(&offset, &n_offset, body, n_body, nwl.distance)) {
        printf("Invalid geometry.\n");
        return -1;
    }

    int n_nodes = (rows + 1) * (cols + 1);
    Node *nodes = malloc(2 * n_nodes * sizeof(Node)); // 2?
    if (!nodes) {
        perror("An error has occurred");
        return -1;
    }

    int node_id = 1;
    for (int i = 0; i <= rows; i++) {
        for (int j = 0; j <= cols; j++) {
            Point p = {X0 + j * cell_size, Y0 + i * cell_size};
            int type = get_point_type(p, offset, n_offset);
            nodes[node_id - 1] = (Node){node_id++, type, p};
        }
    }

    for (int i = 0; i < n_nodes; i++) {
        double frac = 2.01;
        if (nodes[i].type != 0) continue;
        if (is_near_body(&(nodes[i].position), offset, n_offset, cell_size, frac)) {
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
            double frac = 4.0;
            if (is_near_body(&(nodes[i].position), offset, n_offset, cell_size, frac)) {
                nodes[i].type = 2;
            }
        }
    }

    Element *elements = malloc(2 * rows * cols * sizeof(Element)); // 2?
    if (!elements) {
        perror("An error has occurred");
        return -1;
    }

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
                    int n5;
                    pentagons = true;
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

    if (pentagons) split_pentagons(&elements, &n_elements);

    int *offset_nodes;
    int n_offset_nodes = 0;
    get_offset_nodes(&offset_nodes, &n_offset_nodes, nodes, n_nodes, cell_size);
    extrude_near_wall_cells(
        &elements, &n_elements,
        &nodes, &n_nodes,
        body, n_body,
        offset_nodes, n_offset_nodes,
        nwl
    );

    Element *boundaries = malloc(4 * (rows + cols) * sizeof(Element)); // 4?
    if (!boundaries) {
        perror("An error has occurred");
        return -1;
    }

    int n_boundaries = 0;

    // TODO: Define Boundary Elements

    write_mesh_file("mesh.msh", nodes, n_nodes, elements, n_elements, boundaries, n_boundaries);

    free(nodes);
    free(elements);
    free(boundaries);
    free(body);
    free(offset);
    free(offset_nodes);

    printf("Mesh completed. Nodes: %d, Cells: %d\n", n_nodes, n_elements);
    return 0;
}