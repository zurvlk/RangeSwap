#include "graph.h"
#include "bmp.h"

#define _GEN_REV_EDGE_ 1

double energy(Graph *G, int *label, int *I, double T);
double pairwise(double i, double j, double T);
double data(int *I, int a, double i);
int update_labels(Graph *G, int alpha, int *label, int *t, int *I, double T);
void set_capacity(Graph *G, int alpha, int *label, int *I, double T);
void set_all_edge(Graph *G, int height, int width, int label_size, int *I);
void set_single_edge(Graph *G, int height, int width);
