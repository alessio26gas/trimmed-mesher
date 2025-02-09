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

int load_body(const char *filename, Point **body, int *n_points) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("An error has occured while opening the file");
        return -1;
    }

    int count = 0;
    double x, y;
    while (fscanf(file, "%lf,%lf", &x, &y) == 2) {
        count++;
    }

    *body = (Point *)malloc(count * sizeof(Point));
    if (!*body) {
        perror("An error has occurred");
        fclose(file);
        return -1;
    }

    rewind(file);
    int i = 0;
    while (fscanf(file, "%lf,%lf", &(*body)[i].x, &(*body)[i].y) == 2) {
        i++;
    }

    *n_points = count;
    fclose(file);
    return 0;
}

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

    Point *body;
    int n_points;
    if (load_body(argv[1], &body, &n_points) != 0) {
        return -1;
    }

    return 0;
}