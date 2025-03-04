#include "nearwall.h"
#include "point.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846
#define EPSILON 1e-12
#define MAX_ITER 100

double f(double SF, NearWallLayer nwl) {
    if (fabs(SF - 1.0) < EPSILON)
        return nwl.first * nwl.n - nwl.distance;
    return nwl.first * (1 - pow(SF, nwl.n)) / (1 - SF) - nwl.distance;
}

double f_h(double SF, NearWallLayer nwl) {
    return 1 - tanh(SF*(1-1/((double) nwl.n)))/tanh(SF) - nwl.first/nwl.distance;
}

double brent(double (*f)(double, NearWallLayer), NearWallLayer nwl) {
    double a = 1.0;
    double b = 1.0e6;
    double fa = f(a, nwl);
    double fb = f(b, nwl);

    if (fa * fb > 0) {
        printf("\nError: near wall layer settings are not valid.\n");
        exit(1);
    }

    double c = a, fc = fa;
    double d_new = b;
    int iter = 0;

    while (iter < MAX_ITER) {
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        double m = 0.5 * (c - b);

        if (fabs(m) <= EPSILON || fb == 0.0) {
            return b;
        }

        if (fabs(m) > EPSILON && fabs(fa - fb) > EPSILON) {
            double p, q;
            if (a == c) {
                p = 2 * m * fb / (fa - fb);
                q = 1 - fb / (fa - fb);
            } else {
                double r = fb / fc, s = fb / fa, t = fa / fc;
                p = s * (2 * m * t * (t - r) - (b - a) * (r - 1));
                q = (t - 1) * (r - 1) * (s - 1);
            }

            if (p > 0) q = -q;
            p = fabs(p);

            if (2 * p < fmin(3 * m * q - fabs(EPSILON * q), fabs(q * (a - b))))
                d_new = b + p / q;
            else
                d_new = b + m;
        } else {
            d_new = b + m;
        }

        a = b;
        fa = fb;
        b = d_new;
        fb = f(b, nwl);

        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
            c = a;
            fc = fa;
        }

        iter++;
    }

    printf("\nBrent method did not converge after %d iterations.\n", MAX_ITER);
    exit(1);
}

int get_nwl_n(NearWallLayer nwl) {
    if (nwl.last > 0) return 1 + log(nwl.last/nwl.first)/log(nwl.SF);
    if (nwl.SF == 1) return nwl.distance/nwl.first - 1;
    return log(1 - nwl.distance/nwl.first * (1 - nwl.SF))/log(nwl.SF);
}

double get_nwl_distance(NearWallLayer nwl) {
    double sum = 0.0;
    for (int i = 0; i < nwl.n; i++) {
        sum += pow(nwl.SF, i);
    }
    return nwl.first * sum;
}

double get_SF(NearWallLayer nwl) {
    if (nwl.n < 3) return 1.0;
    if (nwl.last > 0) return pow(nwl.last / nwl.first, 1.0 / (nwl.n - 1));
    if (fabs(nwl.distance - nwl.first * nwl.n) < EPSILON) return 1.0;
    return brent(f, nwl);
}

bool compute_offset(Point **offset, int *n_offset, Point *body, int n_body, double nwl_distance, double cell_size) {

    double A = 0.0;
    for (int i = 0; i < n_body; i++) {
        int j = (i + 1) % n_body;
        A += (body[i].x * body[j].y) - (body[j].x * body[i].y);
    }
    int v = 2 * (A < 0) - 1;

    int k = n_body;

    for (int i = 0; i < n_body; i++) {
        int h = (n_body + i - 1) % n_body;
        int j = (i + 1) % n_body;

        double tx_h = (body[h].x - body[i].x);
        double ty_h = (body[h].y - body[i].y);
        double lh = sqrt(tx_h * tx_h + ty_h * ty_h);
        tx_h /= lh;
        ty_h /= lh;

        double tx_j = (body[j].x - body[i].x);
        double ty_j = (body[j].y - body[i].y);
        double lj = sqrt(tx_j * tx_j + ty_j * ty_j);
        tx_j /= lj;
        ty_j /= lj;

        double angle = acos(tx_h * tx_j + ty_h * ty_j);
        double cross_product = tx_h * ty_j - ty_h * tx_j;
        angle += (cross_product * v < 0) * (2*PI - 2*angle);

        if (angle < EPSILON || angle > 2 * PI - EPSILON) return false;
        if (angle < PI) k += 2 * ceil((PI - angle) * nwl_distance / cell_size);
    }

    *n_offset = k;
    *offset = (Point *)malloc((*n_offset) * sizeof(Point));

    k = 0;
    for (int i = 0; i < n_body; i++) {
        int h = (n_body + i - 1) % n_body;
        int j = (i + 1) % n_body;

        double tx_h = (body[h].x - body[i].x);
        double ty_h = (body[h].y - body[i].y);
        double lh = sqrt(tx_h * tx_h + ty_h * ty_h);
        tx_h /= lh;
        ty_h /= lh;

        double tx_j = (body[j].x - body[i].x);
        double ty_j = (body[j].y - body[i].y);
        double lj = sqrt(tx_j * tx_j + ty_j * ty_j);
        tx_j /= lj;
        ty_j /= lj;

        double angle = acos(tx_h * tx_j + ty_h * ty_j);
        double cross_product = tx_h * ty_j - ty_h * tx_j;
        angle += (cross_product * v < 0) * (2*PI - 2*angle);

        if (angle < PI) {
            double nhx = ty_h * v;
            double nhy = -tx_h * v;

            double njx = -ty_j * v;
            double njy = tx_j * v;

            double phih = atan2(nhy, nhx);
            double phij = atan2(njy, njx);

            int N = 1 + 2 * ceil((PI - angle) * nwl_distance / cell_size);

            for (int n = 0; n < N; n++) {
                double alpha = ((double) n) / (N - 1);
                double phin = atan2(
                    (1 - alpha) * sin(phih) + alpha * sin(phij),
                    (1 - alpha) * cos(phih) + alpha * cos(phij)
                );
                (*offset)[k].x = body[i].x + cos(phin) * nwl_distance;
                (*offset)[k].y = body[i].y + sin(phin) * nwl_distance;
                k++;
            }

        } else {
            double bx = tx_h + tx_j;
            double by = ty_h + ty_j;
            double lb = sqrt(bx * bx + by * by);
    
            if (lb > EPSILON) {
                bx /= lb;
                by /= lb;    
            } else {
                bx = -ty_j * v;
                by = tx_j * v;
            }

            (*offset)[k].x = body[i].x + bx * nwl_distance / sin(angle/2);
            (*offset)[k].y = body[i].y + by * nwl_distance / sin(angle/2);
            k++;
        }
    }
    return true;
}

int nearest_node(Node *nodes, int *ids, int n, int current, int *flag, double cell_size) {
    int nearest = -1;
    double min_distance = 2 * cell_size;

    for (int i = 0; i < n; i++) {
        if (!flag[i]) {
            double distance = get_distance(nodes[current].position, nodes[i].position);
            if (distance < min_distance) {
                min_distance = distance;
                nearest = i;
            }
        }
    }
    return nearest;
}

Point nearest_point(Point p, Point *body, int n_points, double nwl_distance) {
    Point nearest;
    double min_distance = 10 * nwl_distance;
    double ABx, ABy, APx, APy, t, Xq, Yq, dx, dy, distance;
    double xp = p.x, yp = p.y;
    int i_min = 0;
    for (int i = 0; i < n_points; i++) {
        int j = (i + 1) % n_points;

        ABx = body[j].x - body[i].x;
        ABy = body[j].y - body[i].y;
        APx = p.x - body[i].x;
        APy = p.y - body[i].y;
        t = (APx * ABx + APy * ABy)/(ABx * ABx + ABy * ABy);
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        Xq = body[i].x + t * ABx;
        Yq = body[i].y + t * ABy;
        dx = Xq - p.x;
        dy = Yq - p.y;
        distance = sqrt(dx * dx + dy * dy);

        if (distance < min_distance) {
            min_distance = distance;
            xp = Xq;
            yp = Yq;
            i_min = i;
        }
    }

    int count = 1;
    for (int i = 0; i < n_points; i++) {
        if (i == i_min) continue;
        int j = (i + 1) % n_points;

        ABx = body[j].x - body[i].x;
        ABy = body[j].y - body[i].y;
        APx = p.x - body[i].x;
        APy = p.y - body[i].y;
        t = (APx * ABx + APy * ABy)/(ABx * ABx + ABy * ABy);
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        Xq = body[i].x + t * ABx;
        Yq = body[i].y + t * ABy;
        dx = Xq - p.x;
        dy = Yq - p.y;
        distance = sqrt(dx * dx + dy * dy);

        if (distance < min_distance + EPSILON) {
            xp += Xq;
            yp += Yq;
            count++;
        }
    }

    nearest.x = xp / count;
    nearest.y = yp / count;
    return nearest;
}

int get_offset_nodes(int **offset_nodes, int *n_offset_nodes, Node *nodes, int n_nodes, double cell_size, double X0, double Y0, int rows, int cols) {

    for (int i = 0; i < n_nodes; i++) {
        if (nodes[i].type == 2) (*n_offset_nodes)++;
    }

    *offset_nodes = malloc(*n_offset_nodes * sizeof(int));
    if (!offset_nodes) {
        perror("An error has occurred");
        return 1;
    }

    Node *type2_nodes = malloc(*n_offset_nodes * sizeof(Node));
    int *ids = malloc(*n_offset_nodes * sizeof(int));
    int count = 0;
    for (int i = 0; i < n_nodes; i++) {
        if (nodes[i].type == 2) {
            type2_nodes[count] = nodes[i];
            ids[count] = nodes[i].id;
            count++;
        }
    }

    int simm = 0;
    for (int i = 0; i < *n_offset_nodes; i++) {
        if (
            type2_nodes[i].position.x < X0 + EPSILON ||
            type2_nodes[i].position.x > X0 + cols * cell_size - EPSILON ||
            type2_nodes[i].position.y < Y0 + EPSILON ||
            type2_nodes[i].position.y > Y0 + rows * cell_size - EPSILON
        ) simm = type2_nodes[i].id;
    }

    int *flag = malloc(*n_offset_nodes * sizeof(int));
    for (int i = 0; i < *n_offset_nodes; i++) flag[i] = 0;

    if (simm != 0) {
        int index = 0;
        int tmp = ids[0];
        Node temp_node = type2_nodes[0];
        for (int i = 0; i < *n_offset_nodes; i++) {
            if (ids[i] == simm) {
                index = i;
                break;
            }
        }
        ids[0] = simm;
        ids[index] = tmp;
        type2_nodes[0] = type2_nodes[index];
        type2_nodes[index] = temp_node;
    }

    int current = 0;
    flag[current] = 1;
    (*offset_nodes)[0] = ids[current];

    for (int i = 1; i < *n_offset_nodes; i++) {
        int nearest = nearest_node(type2_nodes, ids, *n_offset_nodes, current, flag, cell_size);
        if (nearest == -1) break;

        flag[nearest] = 1;
        (*offset_nodes)[i] = ids[nearest];
        current = nearest;
    }

    if (simm != 0) return 1;
    return 0;
}

void extrude_near_wall_cells(
    Element **elements, int *n_elements,
    Node **nodes, int *n_nodes,
    Point *body, int n_body,
    int *offset_nodes, int n_offset_nodes,
    NearWallLayer nwl, int simmetry,
    double X0, int cols, double cell_size
) {
    Point *pA = malloc(n_offset_nodes * sizeof(Point));
    for (int i = 0; i < n_offset_nodes; i++) {
        pA[i] = (*nodes)[offset_nodes[i]-1].position;
    }

    Point *pB = malloc(n_offset_nodes * sizeof(Point));
    for (int i = 0; i < n_offset_nodes; i++) {
        pB[i] = nearest_point(pA[i], body, n_body, nwl.distance);
    }

    int *r = malloc(n_offset_nodes * sizeof(int));
    int *l = malloc(n_offset_nodes * sizeof(int));
    int *c = malloc(n_offset_nodes * sizeof(int));
    for (int i = 0; i < n_offset_nodes; i++) {
        c[i] = 0;
    }
    Point pBold[n_offset_nodes];
    bool singularities = true;

    int n_el = n_offset_nodes;
    int n_el0 = 0;
    if (simmetry == 1) {
        n_el0 = 1;
        n_el = n_offset_nodes - 1;
    }

    int iter = 0;
    while (singularities && iter < nwl.surf_max_iter) {
        iter++;
        singularities = false;
        for (int i = n_el0; i < n_el; i++) {
            int j = (i + 1) % n_offset_nodes;
            int h = (n_offset_nodes + i - 1) % n_offset_nodes;
            r[i] = 0; l[i] = 0;
            if (points_are_close(pB[i], pB[j], nwl.min_surf_distance + EPSILON)) r[i] = 1;
            if (points_are_close(pB[i], pB[h], nwl.min_surf_distance + EPSILON)) l[i] = 1;
        }

        for (int i = n_el0; i < n_el; i++) {
            pBold[i] = pB[i];
        }

        for (int i = n_el0; i < n_el; i++) {
            if (l[i] == 1 && r[i] == 1) {
                singularities = true;
                c[i] = 1;
            }
            if (l[i] == 1 && r[i] == 0) {
                singularities = true;
                int j = (i + 1) % n_offset_nodes;
                if (c[i] == 1) {
                    pB[i].x = 0.5 * pBold[i].x + 0.5 * pBold[j].x;
                    pB[i].y = 0.5 * pBold[i].y + 0.5 * pBold[j].y;                   
                } else {
                    pB[i].x = 2.0 * pBold[i].x / 3.0 + pBold[j].x / 3.0;
                    pB[i].y = 2.0 * pBold[i].y / 3.0 + pBold[j].y / 3.0;
                }
            }
            if (l[i] == 0 && r[i] == 1) {
                singularities = true;
                int h = (n_offset_nodes + i - 1) % n_offset_nodes;
                if (c[i] == 1) {
                    pB[i].x = 0.5 * pBold[i].x + 0.5 * pBold[h].x;
                    pB[i].y = 0.5 * pBold[i].y + 0.5 * pBold[h].y;    
                } else {
                    pB[i].x = 2.0 * pBold[i].x / 3.0 + pBold[h].x / 3.0;
                    pB[i].y = 2.0 * pBold[i].y / 3.0 + pBold[h].y / 3.0;    
                }
            }
        }
    }

    for (int i = 0; i < n_offset_nodes; i++) {
        Point pBtemp;
        do {
            pBtemp = pB[i];
            pB[i] = nearest_point(pB[i], body, n_body, nwl.distance);
        } while (!points_are_equal(pB[i], pBtemp));
    }

    double A = 0.0;
    for (int i = 0; i < n_offset_nodes; i++) {
        int j = (i + 1) % n_offset_nodes;
        A += (pA[i].x * pA[j].y) - (pA[j].x * pA[i].y);
    }
    bool clockwise = (A < 0);

    if (nwl.n == 1) {
        int kk = *n_nodes;

        for (int i = 0; i < n_offset_nodes; i++) {
            (*nodes)[*n_nodes].position = pB[i];
            (*nodes)[*n_nodes].type = 3;
            (*nodes)[*n_nodes].id = *n_nodes + 1;
            (*n_nodes)++;
        }

        for (int i = 0; i < n_el; i++) {
            int n1 = (*nodes)[kk + i].id;
            int n2 = (*nodes)[kk + (i + 1) % n_offset_nodes].id;
            int n3 = offset_nodes[(i + 1) % n_offset_nodes];
            int n4 = offset_nodes[i];

            if (clockwise) {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            (*n_elements)++;
        }
        return;
    }

    if (nwl.n == 2) {
        double x = nwl.first / nwl.distance;
        int kk = *n_nodes;

        for (int i = 0; i < n_offset_nodes; i++) {
            (*nodes)[*n_nodes].position = pB[i];
            (*nodes)[*n_nodes].type = 3;
            (*nodes)[*n_nodes].id = *n_nodes + 1;
            (*n_nodes)++;
        }

        for (int i = 0; i < n_offset_nodes; i++) {
            (*nodes)[*n_nodes].position.x = x*pA[i].x + (1-x)*pB[i].x;
            (*nodes)[*n_nodes].position.y = x*pA[i].y + (1-x)*pB[i].y;
            (*nodes)[*n_nodes].type = 3;
            (*nodes)[*n_nodes].id = *n_nodes + 1;
            (*n_nodes)++;
        }

        for (int i = 0; i < n_el; i++) {
            int n1 = (*nodes)[kk + i].id;
            int n2 = (*nodes)[kk + (i + 1) % n_offset_nodes].id;
            int n3 = (*nodes)[kk + n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n4 = (*nodes)[kk + n_offset_nodes + i].id;

            if (clockwise) {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            (*n_elements)++;
        }

        for (int i = 0; i < n_el; i++) {
            int n1 = (*nodes)[kk + n_offset_nodes + i].id;
            int n2 = (*nodes)[kk + n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n3 = offset_nodes[(i + 1) % n_offset_nodes];
            int n4 = offset_nodes[i];

            if (clockwise) {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            (*n_elements)++;
        }
        return;
    }

    if (nwl.n == 3 && nwl.last > 0) {
        double x[2];
        x[0] = nwl.first / nwl.distance;
        x[1] = 1 - nwl.last / nwl.distance;
        int kk = *n_nodes;

        for (int i = 0; i < n_offset_nodes; i++) {
            (*nodes)[*n_nodes].position = pB[i];
            (*nodes)[*n_nodes].type = 3;
            (*nodes)[*n_nodes].id = *n_nodes + 1;
            (*n_nodes)++;
        }

        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < n_offset_nodes; i++) {
                (*nodes)[*n_nodes].position.x = x[j]*pA[i].x + (1-x[j])*pB[i].x;
                (*nodes)[*n_nodes].position.y = x[j]*pA[i].y + (1-x[j])*pB[i].y;
                (*nodes)[*n_nodes].type = 3;
                (*nodes)[*n_nodes].id = *n_nodes + 1;
                (*n_nodes)++;
            }
        }

        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < n_el; i++) {
                int n1 = (*nodes)[kk + j * n_offset_nodes + i].id;
                int n2 = (*nodes)[kk + j * n_offset_nodes + (i + 1) % n_offset_nodes].id;
                int n3 = (*nodes)[kk + (j + 1) * n_offset_nodes + (i + 1) % n_offset_nodes].id;
                int n4 = (*nodes)[kk + (j + 1) * n_offset_nodes + i].id;

                if (clockwise) {
                    (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n1, n2, n3, n4}};
                } else {
                    (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n4, n3, n2, n1}};
                }
                (*n_elements)++;
            }
        }

        for (int i = 0; i < n_el; i++) {
            int n1 = (*nodes)[kk + n_offset_nodes + i].id;
            int n2 = (*nodes)[kk + n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n3 = offset_nodes[(i + 1) % n_offset_nodes];
            int n4 = offset_nodes[i];

            if (clockwise) {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            (*n_elements)++;
        }
        return;
    }

    double *x = malloc((nwl.n - 1) * sizeof(double));
    double K = 0;
    for (int i = 0; i < nwl.n - 1; i++) {
        K += pow(nwl.SF, i);
        x[i] = nwl.first / nwl.distance * K;
    }

    // HYPERBOLIC TANGENT DISTRIBUTION
    // nwl.SF = brent(f_h, nwl);
    // for (int i = 0; i < nwl.n - 1; i++) {
    //     x[i] = 1 + tanh(nwl.SF * ((i + 1.0)/(nwl.n) - 1))/tanh(nwl.SF);
    // }

    int kk = *n_nodes;

    for (int i = 0; i < n_offset_nodes; i++) {
        (*nodes)[*n_nodes].position = pB[i];
        (*nodes)[*n_nodes].type = 3;
        (*nodes)[*n_nodes].id = *n_nodes + 1;
        (*n_nodes)++;
    }

    for (int j = 0; j < nwl.n - 1; j++) {
        for (int i = 0; i < n_offset_nodes; i++) {
            (*nodes)[*n_nodes].position.x = x[j]*pA[i].x + (1-x[j])*pB[i].x;
            (*nodes)[*n_nodes].position.y = x[j]*pA[i].y + (1-x[j])*pB[i].y;
            (*nodes)[*n_nodes].type = 3;
            (*nodes)[*n_nodes].id = *n_nodes + 1;
            (*n_nodes)++;
        }
    }

    for (int j = 0; j < nwl.n - 1; j++) {
        for (int i = 0; i < n_el; i++) {
            int n1 = (*nodes)[kk + j * n_offset_nodes + i].id;
            int n2 = (*nodes)[kk + j * n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n3 = (*nodes)[kk + (j + 1) * n_offset_nodes + (i + 1) % n_offset_nodes].id;
            int n4 = (*nodes)[kk + (j + 1) * n_offset_nodes + i].id;

            if (clockwise) {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n1, n2, n3, n4}};
            } else {
                (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n4, n3, n2, n1}};
            }
            (*n_elements)++;
        }
    }

    for (int i = 0; i < n_el; i++) {
        int n1 = (*nodes)[kk + (nwl.n - 1) * n_offset_nodes + i].id;
        int n2 = (*nodes)[kk + (nwl.n - 1) * n_offset_nodes + (i + 1) % n_offset_nodes].id;
        int n3 = offset_nodes[(i + 1) % n_offset_nodes];
        int n4 = offset_nodes[i];

        if (clockwise) {
            (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n1, n2, n3, n4}};
        } else {
            (*elements)[*n_elements] = (Element){*n_elements + 1, 3, 4, {n4, n3, n2, n1}};
        }
        (*n_elements)++;
    }
}