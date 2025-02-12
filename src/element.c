#include "element.h"
#include <stdio.h>
#include <stdlib.h>

void tria_vertices(int *n1, int *n2, int *n3, int *n4, Point *body, int n_points, Node **nodes, int *n_nodes) {
    if ((*nodes)[*n1 - 1].type == 0) {
        *n3 = *n4;
        if ((*nodes)[*n2 - 1].type != 2)
        find_and_update(n2, n1, body, n_points, nodes, n_nodes);
        if ((*nodes)[*n3 - 1].type != 2)
        find_and_update(n3, n1, body, n_points, nodes, n_nodes);
    } else if ((*nodes)[*n2 - 1].type == 0) {
        if ((*nodes)[*n3 - 1].type != 2)
        find_and_update(n3, n2, body, n_points, nodes, n_nodes);
        if ((*nodes)[*n1 - 1].type != 2)
        find_and_update(n1, n2, body, n_points, nodes, n_nodes);
    } else if ((*nodes)[*n3 - 1].type == 0) {
        *n1 = *n2;
        *n2 = *n3;
        *n3 = *n4;
        if ((*nodes)[*n1 - 1].type != 2)
        find_and_update(n1, n2, body, n_points, nodes, n_nodes);
        if ((*nodes)[*n3 - 1].type != 2)
        find_and_update(n3, n2, body, n_points, nodes, n_nodes);
    } else if ((*nodes)[*n4 - 1].type == 0) {
        *n2 = *n3;
        *n3 = *n4;
        if ((*nodes)[*n1 - 1].type != 2)
        find_and_update(n1, n3, body, n_points, nodes, n_nodes);
        if ((*nodes)[*n2 - 1].type != 2)
        find_and_update(n2, n3, body, n_points, nodes, n_nodes);
    }
}

void quad_vertices(int *n1, int *n2, int *n3, int *n4, Point *body, int n_points, Node **nodes, int *n_nodes) {
    if ((*nodes)[*n1 - 1].type != 0 && (*nodes)[*n2 - 1].type != 0) {
        if ((*nodes)[*n1 - 1].type != 2)
        find_and_update(n1, n4, body, n_points, nodes, n_nodes);
        if ((*nodes)[*n2 - 1].type != 2)
        find_and_update(n2, n3, body, n_points, nodes, n_nodes);
    } else if ((*nodes)[*n2 - 1].type != 0 && (*nodes)[*n3 - 1].type != 0) {
        if ((*nodes)[*n2 - 1].type != 2)
        find_and_update(n2, n1, body, n_points, nodes, n_nodes);
        if ((*nodes)[*n3 - 1].type != 2)
        find_and_update(n3, n4, body, n_points, nodes, n_nodes);
    } else if ((*nodes)[*n3 - 1].type != 0 && (*nodes)[*n4 - 1].type != 0) {
        if ((*nodes)[*n3 - 1].type != 2)
        find_and_update(n3, n2, body, n_points, nodes, n_nodes);
        if ((*nodes)[*n4 - 1].type != 2)
        find_and_update(n4, n1, body, n_points, nodes, n_nodes);
    } else if ((*nodes)[*n4 - 1].type != 0 && (*nodes)[*n1 - 1].type != 0) {
        if ((*nodes)[*n4 - 1].type != 2)
        find_and_update(n4, n3, body, n_points, nodes, n_nodes);
        if ((*nodes)[*n1 - 1].type != 2)
        find_and_update(n1, n2, body, n_points, nodes, n_nodes);
    }
}

int penta_vertices(int *n1, int *n2, int *n3, int *n4, int *n5, Point *body, int n_points, Node **nodes, int *n_nodes) {
    if ((*nodes)[*n1 - 1].type == 1) {
        *n5 = *n4;
        *n4 = *n3;
        *n3 = *n2;
        find_and_update(n2, n1, body, n_points, nodes, n_nodes);
        find_and_update(n1, n5, body, n_points, nodes, n_nodes);
        return 1;
    } else if ((*nodes)[*n2 - 1].type == 1) {
        *n5 = *n4;
        *n4 = *n3;
        find_and_update(n3, n2, body, n_points, nodes, n_nodes);
        find_and_update(n2, n1, body, n_points, nodes, n_nodes);
        return 2;
    } else if ((*nodes)[*n3 - 1].type == 1) {
        *n5 = *n4;
        find_and_update(n4, n3, body, n_points, nodes, n_nodes);
        find_and_update(n3, n2, body, n_points, nodes, n_nodes);
        return 3;
    } else if ((*nodes)[*n4 - 1].type == 1) {
        *n5 = *n1;
        find_and_update(n5, n4, body, n_points, nodes, n_nodes);
        find_and_update(n4, n3, body, n_points, nodes, n_nodes);
        return 4;
    }
}

void split_pentagons(Element **elements, int *num_elements) {
    int new_count = *num_elements;

    for (int i = 0; i < *num_elements; i++) {
        if ((*elements)[i].type == 4) new_count += 2;
    }

    Element *new_elements = (Element *)malloc(2 * new_count * sizeof(Element));
    int new_id = 1, new_index = 0;

    for (int i = 0; i < *num_elements; i++) {
        Element e = (*elements)[i];

        if (e.type == 4) {
            int ci = (e.flag + 2) % 5;

            for (int j = 0; j < 3; j++, new_index++) {
                new_elements[new_index] = (Element){ 
                    .id = new_id++, .type = 2, .n_nodes = 3, 
                    .nodes = { e.nodes[ci], e.nodes[(ci + j + 1) % 5], e.nodes[(ci + j + 2) % 5] }
                };
            }
        } else {
            new_elements[new_index] = e;
            new_elements[new_index].id = new_id++;
            new_index++;
        }
    }

    free(*elements);
    *elements = new_elements;
    *num_elements = new_count;
}