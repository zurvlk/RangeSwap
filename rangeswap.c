#include "bmp.h"
#include "rangeswap.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#define INF DBL_MAX

int isin_array(int *array, int target, int size) {
    int high, mid, low;
    
    high = size;
    low = 1;

    // while (low <= high) {
    //     mid = (low + high) / 2;

    //     /* 値が見つかればループを抜ける */
    //     if (array[mid] == target) {
    //       return 1;
    //       /* 値の大小を調べて探索範囲を狭める */
    //     } else if (array[mid] < target) {
    //       low = mid + 1;
    //     } else {
    //       high = mid - 1;
    //     }
    // }

    for (int i = 1; i <= size; i++) {
        if(array[i] == target) return 1;
    }
    return 0;
}

// 2つの配列の値がすべて同一 0
// 2つの配列の値が異なる 1
int cmparray(int *array1, int *array2, int size) {
    int i;

    for (i = 1; i <= size; i++) {
        if(array1[i] != array2[i]) return 1;
    }
    return 0;
}

void cpyarray(int *terget, int *source, int size) {
    int i;

    for (i = 1; i <= size; i++) {
        terget[i] = source[i];
    }
}

double h(double n) {
    /*
    return n > 0 ? n : -n;
    /*/
    return n * n;
    // */
}

double dabs(double a, double b) {
    return a - b > 0 ? a - b : b - a;
}
double fmin(double i, double j) {
    return i < j ? i : j;
}

double pairwise(double i, double j, double T) {
    /*
    return fmin(1.0 * h(i - j), T);
    /*/
    return h(i - j);
    // */
}

double data(int *I, int a, double i) {
    return 1.0 * dabs(I[a], i);
}


double energy(Graph *G, int *label, int *I, double T) {
    int i;
    double energy = 0;
    //* Dterm
    for (i = 1; i <= G->n - 2; i++) {
        energy += data(I, i, label[i]);
    }
    // */
    // Vterm
    for (i = 1; i <= G->m - 2 * (G->n - 2); i++) {
        energy += pairwise(label[G->tail[i]], label[G->head[i]], T);
    }
    return energy;
}

double r(int i, int k, int grids_edge) {
    return grids_edge * (- 0.5) * (h(k + 1 - i) + h(i + 1));
}

double e_cost(int i, int j) {
    double cost = 0.5 * (h(i - j + 1)- 2 * h(i - j) + h(i - j - 1));
    return cost;
}

void set_all_edge(Graph *G, int height, int width, int label_size, int *I) {
    int i, j, k, l;
    int tail, head, t_base, h_base, grids_node, grids_edge, source, sink, edge_count, current_edge;
    int s2i_begin, i2t_begin, depth_begin;
    double *min, r0, rk;



    if (((min = (double *) malloc(sizeof(double) * G->n))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [min = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 0; i < G->n; i++) min[i] = INF;

    // 格子部分1階層分の点数合計
    grids_node = height * width;

    // 格子部分1階層分の枝数合計
    grids_edge = (height - 1) * width + height * (width - 1);

    for (i = 1; i < G->n; i++) G->capa[i] = 0;
    source = grids_node * label_size + 1;
    sink = source + 1;


    setSource(G, source);
    setSink(G, sink);

    edge_count = 1;
    // source->i1
    s2i_begin = edge_count;
    r0 = r(0, label_size, grids_edge);
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, source, i, 0);
        G->capa[edge_count] = data(I, i, 0) - r0 + e_cost(0, 1);
        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }
    
    // source->i2~ik
    for (i = grids_node + 1; i <= grids_node * label_size; i++) {
        setEdge(G, edge_count, source, i, 0);
        G->capa[edge_count] = e_cost(0, (i - 1) / grids_node + 1);
        edge_count++;
    }

    // i1~ik-1->sink
    for (i = 1; i <= (label_size - 1) * grids_node; i++) {
        setEdge(G, edge_count, i, sink, 0);
        G->capa[edge_count] = e_cost((i - 1) / grids_node + 1, label_size + 1);
        edge_count++;
    }
    // ik->sink
    i2t_begin = edge_count;
    rk = r(label_size , label_size, grids_edge);
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, i + grids_node * (label_size - 1), sink, 0);
        G->capa[edge_count] = data(I, i, label_size) - rk + e_cost(label_size, label_size + 1);
        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }

    // depth
    depth_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < label_size; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, tail, head, 0);
            G->capa[edge_count] = data(I, i, j) - r(j, label_size, grids_edge);
            if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
            edge_count++;
            tail = head;
        }
    }

    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 0; k < label_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < label_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    G->capa[edge_count] = e_cost(k + 1, l + 1);
                    edge_count++;
                }
            }
        }
    }
    
    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 0; k < label_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < label_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    G->capa[edge_count] = e_cost(k + 1, l + 1);
                    edge_count++;
                }
            }
        }
    }


    // reverce edge
    // depth
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < label_size; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, head, tail, INF);
            edge_count++;
            tail = head;
        }
    }

    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 0; k < label_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < label_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    G->capa[edge_count] = e_cost(k + 1, l + 1);
                    edge_count++;
                }
            }
        }
    }

    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 0; k < label_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < label_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    G->capa[edge_count] = e_cost(k + 1, l + 1);
                    edge_count++;
                }
            }
        }
    }


    // each node 2 source 
    for (i = 1; i <= grids_node * label_size; i++) {
        setEdge(G, edge_count, i, source, 0);
        G->capa[edge_count] = e_cost(0, (i - 1) / grids_node + 1);
        edge_count++;
    }
     // each node 2 sink
    for (i = 1; i <= grids_node * label_size; i++) {
        setEdge(G, edge_count, sink, i, 0);
        G->capa[edge_count] = e_cost((i - 1) / grids_node + 1, label_size + 1);
        edge_count++;
    }


    //  s->tの一連の枝から定数値を引く処理
    for (i = s2i_begin; i <= grids_node; i++) {
        G->capa[i] -= min[i];
    }
    for (i = i2t_begin; i <= grids_node; i++) {
        G->capa[i + i2t_begin - 1] -= min[i];
    }
    current_edge = depth_begin;
    for (i = 1; i <= grids_node; i++) {
        for (j = 1; j < label_size; j++) {
            G->capa[current_edge] -=min[i];
            current_edge++;
        }
    }

    free(min);
    // printf("total edge : %d\n", edge_count - 1);
    return;
}

//枝を作る関数makeedge(グラフ,高さ,幅)
void set_single_edges(Graph *G, int height, int width) {
    int i, j, edge_count;
    int tail, head, source, sink;

    source = G->n - 1;
    sink = source + 1;
    setSource(G, source);
    setSink(G, sink);
    
    edge_count = 1;
    //点と点の間の枝（横）
    for (i = 1; i < height + 1; i++) {
        for (j = 1; j < width; j++) {
            tail =  (i - 1) * width + j;
            head =  tail + 1;
            setEdge(G, edge_count, tail, head, 0);
            edge_count++;
        }
    }
 
    //点と点の間の枝（縦）
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            tail = (i - 1) * width + j;
            head = tail + width;
            setEdge(G, edge_count, tail, head, 0);
            edge_count++;
        }
    }
 
    //sourceと点の間の枝
    for (i = 1; i < height * width + 1; i++){
       setEdge(G, edge_count, G->src, i, 0);
       edge_count++;
    }
    
    //点とsinkの間の枝
    for (i = 1; i < height * width + 1; i++){
       setEdge(G, edge_count, i, G->sink, 0);
       edge_count++;
    }
    return;
 }
int make_label_index(Graph *G, int *label, int *label_index, int alpha, int beta) {
    int i, arraysize;

    arraysize = 1;
    for (i = 1; i <= G->n - 2; i++) {
        if (label[i] <= beta && label[i] >= alpha) {
            label_index[arraysize] = i;
            arraysize++;
        }
    }
    return arraysize;
}
 // set_edge for rangeswap
void set_edge(Graph *G, int height, int width, int alpha, int beta, int *label_index, int size, int *I) {
    int i, j, k, l;
    int tail, head, t_base, h_base, grids_node, grids_edge, source, sink, edge_count, current_edge;
    int s2i_begin, i2t_begin, depth_begin, label_size;
    double *min, r0, rk;


    label_size = beta - alpha;
    if (((min = (double *) malloc(sizeof(double) * G->n))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [min = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 0; i < G->n; i++) min[i] = INF;

    // 格子部分1階層分の点数合計
    grids_node = height * width;

    // 格子部分1階層分の枝数合計
    grids_edge = (height - 1) * width + height * (width - 1);

    for (i = 1; i < G->n; i++) G->capa[i] = 0;
    source = grids_node * label_size + 1;
    sink = source + 1;


    setSource(G, source);
    setSink(G, sink);

    edge_count = 1;
    // source->i1
    s2i_begin = edge_count;
    r0 = r(0, label_size, grids_edge);
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, source, i, 0);

        if(isin_array(label_index, i, size)) {
            G->capa[edge_count] = data(I, i, 0) - r0 + e_cost(0, 1);
        }
        
        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }
    
    // source->i2~ik
    for (i = grids_node + 1; i <= grids_node * label_size; i++) {
        setEdge(G, edge_count, source, i, 0);

        if(isin_array(label_index, i % grids_node, size)) {
            G->capa[edge_count] = e_cost(0, (i - 1) / grids_node + 1);
        }

        edge_count++;
    }

    // i1~ik-1->sink
    for (i = 1; i <= (label_size - 1) * grids_node; i++) {
        setEdge(G, edge_count, i, sink, 0);

        if(isin_array(label_index, i % grids_node, size)) {
            G->capa[edge_count] = e_cost((i - 1) / grids_node + 1, label_size + 1);
        }

        edge_count++;
    }
    // ik->sink
    i2t_begin = edge_count;
    rk = r(label_size , label_size, grids_edge);
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, i + grids_node * (label_size - 1), sink, 0);

        if(isin_array(label_index, i % grids_node, size)) {
            G->capa[edge_count] = data(I, i, label_size) - rk + e_cost(label_size, label_size + 1);
        }

        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }

    // depth
    depth_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < label_size; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, tail, head, 0);

            if(isin_array(label_index, i , size)) {
                G->capa[edge_count] = data(I, i, j) - r(j, label_size, grids_edge);
            }

            if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
            edge_count++;
            tail = head;
        }
    }

    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 0; k < label_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < label_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        G->capa[edge_count] = e_cost(k + 1, l + 1);
                    }
                    edge_count++;
                }
            }
        }
    }
    
    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 0; k < label_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < label_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        G->capa[edge_count] = e_cost(k + 1, l + 1);
                    }
                    edge_count++;
                }
            }
        }
    }


    // reverce edge
    // depth
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < label_size; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, head, tail, 0);
            if(isin_array(label_index, i, size)) {
                G->capa[edge_count] = INF;
            }
            edge_count++;
            tail = head;
        }
    }

    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 0; k < label_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < label_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        G->capa[edge_count] = e_cost(k + 1, l + 1);
                    }
                    edge_count++;
                }
            }
        }
    }

    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 0; k < label_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < label_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        G->capa[edge_count] = e_cost(k + 1, l + 1);
                    }
                    edge_count++;
                }
            }
        }
    }


    // each node 2 source 
    for (i = 1; i <= grids_node * label_size; i++) {
        setEdge(G, edge_count, i, source, 0);
        if(isin_array(label_index, i % grids_node, size)) {
            G->capa[edge_count] = e_cost(0, (i - 1) / grids_node + 1);
        }
        edge_count++;
    }
     // each node 2 sink
    for (i = 1; i <= grids_node * label_size; i++) {
        setEdge(G, edge_count, sink, i, 0);
        if(isin_array(label_index, i % grids_node, size)) {
            G->capa[edge_count] = e_cost((i - 1) / grids_node + 1, label_size + 1);
        }
        edge_count++;
    }


    //  s->tの一連の枝から定数値を引く処理
    for (i = s2i_begin; i <= grids_node; i++) {
        G->capa[i] -= min[i];
    }
    for (i = i2t_begin; i <= grids_node; i++) {
        G->capa[i + i2t_begin - 1] -= min[i];
    }
    current_edge = depth_begin;
    for (i = 1; i <= grids_node; i++) {
        for (j = 1; j < label_size; j++) {
            G->capa[current_edge] -=min[i];
            current_edge++;
        }
    }

    free(min);
    // printf("total edge : %d\n", edge_count - 1);
    return;
}