#include <stdio.h>
#include <stdlib.h>

#include "main.h"
#include "input.h"
#include "export.h"

void write_mesh_file() {
    FILE *file = fopen(input.outputfile, "w");
    if (!file) {
        perror("Errore nell'apertura del file di output");
        return;
    }

    fprintf(file, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    fprintf(file, "$Nodes\n");
    fprintf(file, "%d\n", n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        fprintf(file, "%d %.12lf %.12lf 0.0\n", nodes[i].id, nodes[i].position.x, nodes[i].position.y);
    }
    fprintf(file, "$EndNodes\n");

    fprintf(file, "$Elements\n");
    fprintf(file, "%d\n", n_elements + n_boundaries);
    for (int i = 0; i < n_elements; i++) {
        if (elements[i].type == 0) continue;
        fprintf(file, "%d %d 2 99 0", elements[i].id, elements[i].type);
        for (int j = 0; j < elements[i].n_nodes; j++) {
            fprintf(file, " %d", elements[i].nodes[j]);
        }
        fprintf(file, "\n");
    }

    for (int i = 0; i < n_boundaries; i++) {
        if (boundaries[i].type == 0) continue;
        fprintf(file, "%d %d 2 99 %d", boundaries[i].id, boundaries[i].type, boundaries[i].flag);
        for (int j = 0; j < boundaries[i].n_nodes; j++) {
            fprintf(file, " %d", boundaries[i].nodes[j]);
        }
        fprintf(file, "\n");
    }

    fprintf(file, "$EndElements\n");

    fclose(file);
}