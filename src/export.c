#include "export.h"
#include <stdio.h>
#include <stdlib.h>

void write_mesh_file(const char *filename, Node *nodes, int n_nodes, Element *elements, int n_elements, Element *boundaries, int n_boundaries) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Errore nell'apertura del file di output");
        return;
    }

    fprintf(file, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    // TODO: Physical Entities

    fprintf(file, "$Nodes\n");
    fprintf(file, "%d\n", n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        fprintf(file, "%d %lf %lf 0.0\n", nodes[i].id, nodes[i].position.x, nodes[i].position.y);
    }
    fprintf(file, "$EndNodes\n");

    fprintf(file, "$Elements\n");
    fprintf(file, "%d\n", n_elements);
    for (int i = 0; i < n_elements; i++) {
        if (elements[i].type == 0) continue;
        fprintf(file, "%d %d 1 99", elements[i].id, elements[i].type);
        for (int j = 0; j < elements[i].n_nodes; j++) {
            fprintf(file, " %d", elements[i].nodes[j]);
        }
        fprintf(file, "\n");
    }

    // TODO: Write Boundary Elements

    fprintf(file, "$EndElements\n");

    fclose(file);
}