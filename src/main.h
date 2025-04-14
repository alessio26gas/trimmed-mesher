#ifndef MAIN_H
#define MAIN_H

#include "node.h"
#include "element.h"

extern Point *body;
extern int n_body;

extern Point *offset;
extern int n_offset;

extern Node *nodes;
extern int n_nodes;

extern Element *elements;
extern int n_elements;

extern Element *boundaries;
extern int n_boundaries;

#endif