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
    int type; // 0 = not enclosed, 1 = enclosed, 2 = near-body
    Point position;
} Node;

typedef struct {
    int id;
    int type;  // 2 = triangle, 3 = quadrilateral, 4 = pentagon
    int n_nodes; // 3 = triangle, 4 = quadrilateral, 5 = pentagon
    int nodes[5];
    int flag; // n1=1 n2=2 n3=3 n4=4
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

int is_enclosed(Point p, Point *body, int n_points) {
    int inside = 0;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;
        if (((body[i].y > p.y) != (body[j].y > p.y)) &&
            (p.x < (body[j].x - body[i].x) * (p.y - body[i].y) / (body[j].y - body[i].y) + body[i].x)) {
            inside = !inside;
        }
    }
    return inside;
}

bool is_near_body(Point *p, Point *body, int n_points, double cell_size) {
    double min_distance = cell_size;
    double ABx, ABy, APx, APy, ABAB, ABAP, t, Xq, Yq, dx, dy, distance;
    double xp = (*p).x, yp = (*p).y;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;

        ABx = body[j].x - body[i].x;
        ABy = body[j].y - body[i].y;
        APx = (*p).x - body[i].x;
        APy = (*p).y - body[i].y;
        t = (APx * ABx + APy * ABy)/(ABx * ABx + ABy * ABy);
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        Xq = body[i].x + t * ABx;
        Yq = body[i].y + t * ABy;
        dx = Xq - (*p).x;
        dy = Yq - (*p).y;
        distance = sqrt(dx * dx + dy * dy);

        if (distance < min_distance) {
            min_distance = distance;
            xp = Xq;
            yp = Yq;
        }
    }
    if (min_distance < cell_size / 4) {
        (*p).x = xp;
        (*p).y = yp;
        return true;
    }
    return false;
}

bool get_intersection(Point a, Point b, Point c, Point d, Point *intersect) {
    double det = (b.x - a.x) * (d.y - c.y) - (b.y - a.y) * (d.x - c.x);
    if (fabs(det) < EPSILON) return false;
    double t = ((c.x - a.x) * (d.y - c.y) - (c.y - a.y) * (d.x - c.x)) / det;
    double u = ((c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x)) / det;
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        intersect->x = a.x + t * (b.x - a.x);
        intersect->y = a.y + t * (b.y - a.y);
        return true;
    }
    return false;
}

bool points_are_equal(Point a, Point b) {
    return (fabs(a.x - b.x) < EPSILON) && (fabs(a.y - b.y) < EPSILON);
}

int find_or_add_node(Point p, Node *nodes, int *n_nodes) {
    for (int i = 0; i < *n_nodes; i++) {
        if (points_are_equal(nodes[i].position, p)) {
            return nodes[i].id;
        }
    }
    nodes[*n_nodes] = (Node){*n_nodes + 1, 0, p};
    return ++(*n_nodes);
}

void find_and_update(int *nA, int *nB, Point *body, int n_points, Node **nodes, int *n_nodes) {
    Point p;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;
        if (get_intersection((*nodes)[*nA-1].position, (*nodes)[*nB-1].position, body[i], body[j], &p)) {
            break;
        }
    }
    *nA = find_or_add_node(p, *nodes, n_nodes);
}

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

    Element *new_elements = (Element *)malloc(new_count * sizeof(Element));
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

void write_mesh_file(const char *filename, Node *nodes, int n_nodes, Element *elements, int n_elements) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Errore nell'apertura del file di output");
        return;
    }

    fprintf(file, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

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

    fprintf(file, "$EndElements\n");

    fclose(file);
}

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

    double X0, Y0;
    if (argc == 7) {
        X0 = atof(argv[5]);
        Y0 = atof(argv[6]);
    } else if (argc == 6) {
        X0 = atof(argv[5]);
        Y0 = - cell_size * rows / 2;
    } else {
        X0 = - cell_size * cols / 2;
        Y0 = - cell_size * rows / 2;
    }

    Point *body;
    int n_points;
    if (load_body(argv[1], &body, &n_points) != 0) {
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
            int type = is_enclosed(p, body, n_points);
            nodes[node_id - 1] = (Node){node_id++, type, p};
        }
    }

    for (int i = 0; i < n_nodes; i++) {
        if (is_near_body(&(nodes[i].position), body, n_points, cell_size)) {
            nodes[i].type = 2;
        }
    }

    Element *elements = malloc(2 * rows * cols * sizeof(Element)); // 2?
    if (!elements) {
        perror("An error has occurred");
        return -1;
    }

    int element_id = 1, n_elements = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int n1 = i * (cols + 1) + j + 1;
            int n2 = n1 + 1;
            int n3 = n1 + (cols + 1) + 1;
            int n4 = n1 + (cols + 1);

            int ex = (nodes[n1-1].type==0) + (nodes[n2-1].type==0) + (nodes[n3-1].type==0) + (nodes[n4-1].type==0);
            int in = (nodes[n1-1].type==1) + (nodes[n2-1].type==1) + (nodes[n3-1].type==1) + (nodes[n4-1].type==1);

            if (in == 0) {
                elements[n_elements++] = (Element){element_id++, 3, 4, {n1, n2, n3, n4}};
                continue;
            }

            switch (ex) {
                case 3:
                    int n5;
                    int flag = penta_vertices(&n1, &n2, &n3, &n4, &n5, body, n_points, &nodes, &n_nodes);
                    elements[n_elements++] = (Element){element_id++, 4, 5, {n1, n2, n3, n4, n5}, flag};
                    break;
                case 2:
                    quad_vertices(&n1, &n2, &n3, &n4, body, n_points, &nodes, &n_nodes);
                    elements[n_elements++] = (Element){element_id++, 3, 4, {n1, n2, n3, n4}};
                    break;
                case 1:
                    tria_vertices(&n1, &n2, &n3, &n4, body, n_points, &nodes, &n_nodes);
                    elements[n_elements++] = (Element){element_id++, 2, 3, {n1, n2, n3}};
                    break;
            }
        }
    }

    split_pentagons(&elements, &n_elements);

    write_mesh_file("mesh.msh", nodes, n_nodes, elements, n_elements);

    free(nodes);
    free(elements);
    free(body);

    printf("Mesh completed. Nodes: %d, Cells: %d\n", n_nodes, n_elements);
    return 0;
}