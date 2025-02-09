#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define EPSILON 1e-8

typedef struct {
    double x, y;
} Point;

typedef struct {
    int id;
    int type;
    Point position;
} Node;

typedef struct {
    int id;
    int type;
    int num_nodes;
    Node nodes[5];
    int pentaflag;
} Element;

int main(int argc, char *argv[]) {

    if (argc < 4 || argc > 5) {
        printf("Usage: %s <body.csv> <cellsize> <rows> [columns]\n", argv[0]);
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
        printf("Usage: %s <body.csv> <cellsize> <rows> [columns]\n", argv[0]);
        return -1;
    }

    return 0;
}