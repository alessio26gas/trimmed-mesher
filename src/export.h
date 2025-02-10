#ifndef EXPORT_H
#define EXPORT_H

#include "node.h"
#include "element.h"

void write_mesh_file(const char *filename, Node *nodes, int n_nodes, Element *elements, int n_elements, Element *boundaries, int n_boundaries);

#endif