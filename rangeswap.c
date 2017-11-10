#include "bmp.h"
#include "rangeswap.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#define INF DBL_MAX

int isin_array(int *array, int target, int size) {
    for (int i = 1; i <= size; i++) {
        if(array[i] > target) break;
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

    if(!function) return n > 0 ? n : -n;
    else return n * n;

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

double d_p(int label, int i, int alpha) {
    return dabs(label, (i + alpha));
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
    for (i = 1; i <= G->n - 2; i++) label_index[i] = 0;
    arraysize = 1;
    for (i = 1; i <= G->n - 2; i++) {
        if (label[i] <= beta && label[i] >= alpha) {
            label_index[arraysize] = i;
            arraysize++;
        }
    }
    return arraysize;
}

int r_4_rangeswap(int i, int li, int height, int width, int label_size, int *label_index, int *label, int index_size, int alpha, int beta) {
    int grids_node = height * width;

    double r_total = 0;
    double k = dabs(alpha, beta);

    if (i > width) {
        // 画素が一番上の行に存在しないとき(iの上が空白でないとき)
        if (isin_array(label_index, i - width, index_size)){
            // iの上の点がalpha-beta間に含まれる
            r_total -= r(li, k, 1);
        } else {
            // iの上の点がalpha-beta間に含まれない-= r(li, k, 1);
            r_total += h(li + alpha - label[i - width]);
        }
    }
    

    if (i < grids_node - width) {
        // 画素が一番下の行に存在しないとき(iの下が空白でないとき)
        if (isin_array(label_index, i + width, index_size)){
            // iの下の点がalpha-beta間に含まない
            r_total -= r(li, k, 1);
        } else {
            // iの下の点がalpha-beta間に含まれない
            r_total += h(li + alpha -  label[i + width]);
        }
    }

    if ((i % width) != 1) {
        // 画素が一番左の列に存在しないとき(iの左が空白でないとき)
        if (isin_array(label_index, i - 1, index_size)){
            // iの左の点がalpha-beta間に含まれる
            r_total -= r(li, k, 1);
        } else {
            // iの左の点がalpha-beta間に含まれない
            r_total += h(li + alpha -  label[i - 1]);
        }
    }

    if ((i % width) != 0) {
        // 画素が一番右の列に存在しないとき(iの右が空白でないとき)
        if (isin_array(label_index, i + 1, index_size)){
            // iの右の点がalpha-beta間に含まれる
            r_total -= r(li, k, 1);
        } else {
            // iの右の点がalpha-beta間に含まれない
            r_total += h(li + alpha - label[i + 1]);
        }
    }
    return r_total;
}

// set_edge for rangeswap
void set_edge(Graph *G, int height, int width, int alpha, int beta, int label_size, int *label, int *label_index, int size, int *I) {
    int i, j, k, l, node;
    int tail, head, t_base, h_base, grids_node, source, sink, edge_count, current_edge;
    int s2i_begin, i2t_begin, depth_begin, range_size, sink_label;
    double *min, r_total;


    if (((min = (double *) malloc(sizeof(double) * G->n))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [min = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 0; i < G->n; i++) min[i] = INF;

    // 格子部分1階層分の点数合計
    grids_node = height * width;
    range_size = beta - alpha;

    for (i = 1; i < G->n; i++) G->capa[i] = 0;
    source = grids_node * range_size + 1;
    sink = source + 1;
    sink_label = range_size + 1;

    setSource(G, source);
    setSink(G, sink);

    edge_count = 1;
    // source->i1
    s2i_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, source, i, 0);
        if(isin_array(label_index, i, size)) {
            r_total = r_4_rangeswap(i, 0, height, width, label_size, label_index, label, size, alpha, beta);
            // G->capa[edge_count] = data(I, i, alpha) - r0 + e_cost(alpha, alpha + 1);
            G->capa[edge_count] = d_p(label[i], 0, alpha) + r_total + e_cost(0, 1);
        }
        
        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }
    // printf("%lf\n", G->capa[edge_count - 1]);
    // source->i2~ik
    for (i = grids_node + 1; i <= grids_node * range_size; i++) {
        setEdge(G, edge_count, source, i, 0);
        node = i % grids_node == 0 ? grids_node : i % grids_node;
        if(isin_array(label_index, node, size)) {
            G->capa[edge_count] = e_cost(0, (i - 1) / grids_node + 1);
        }

        edge_count++;
    }
    // printf("%d %d \n",i,  i % grids_node == 0 ? grids_node : i % grids_node);
    // printf("::%d\n", range_size); 

    // i1~ik-1->sink
    for (i = 1; i <= (range_size - 1) * grids_node; i++) {
        setEdge(G, edge_count, i, sink, 0);
        node = i % grids_node == 0 ? i : i % grids_node;
        if(isin_array(label_index, node, size)) {
            G->capa[edge_count] = e_cost((i - 1) / grids_node + 1, sink_label);
        }

        edge_count++;
    }
    // printf("%lf\n", G->capa[edge_count - 1]);
    // ik->sink
    i2t_begin = edge_count;
    // rk = r(beta , label_size, grids_edge);
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, i + grids_node * (range_size - 1), sink, 0);
        node = i % grids_node == 0 ? grids_node : i % grids_node;
        if(isin_array(label_index, node, size)) {
            r_total = r_4_rangeswap(i, beta, height, width, label_size, label_index, label, size, alpha, beta);
            G->capa[edge_count] = d_p(I[i], range_size , alpha) + r_total + e_cost(range_size, sink_label);
        }

        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }

    // depth
    depth_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < range_size; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, tail, head, 0);

            if(isin_array(label_index, i , size)) {
                r_total = r_4_rangeswap(i, j, height, width, label_size, label_index, label, size, alpha, beta);
                G->capa[edge_count] = d_p(I[i], j, alpha) + r_total;
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
            for (k = 0; k < range_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < range_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        // head, tail in label_index
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
            for (k = 0; k < range_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < range_size; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(label_index, t_base, size) && isin_array(label_index, h_base, size)) {
                        G->capa[edge_count] = e_cost(k + 1, l + 1);
                        // head, tail in label_index
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
        for (j = 1; j < range_size; j++) {
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
            for (k = 0; k < range_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < range_size; l++) {
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
            for (k = 0; k < range_size; k++) {
                tail = t_base + k * grids_node;
                for (l = 0; l < range_size; l++) {
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
    for (i = 1; i <= grids_node * range_size; i++) {
        setEdge(G, edge_count, i, source, 0);
        node = i % grids_node == 0 ? i : i % grids_node;
        if(isin_array(label_index, node, size)) {
            G->capa[edge_count] = e_cost(0, (i - 1) / grids_node + 1);
        }
        edge_count++;
    }
     // each node 2 sink
    for (i = 1; i <= grids_node * range_size; i++) {
        setEdge(G, edge_count, sink, i, 0);
        node = i % grids_node == 0 ? i : i % grids_node;
        if(isin_array(label_index, node, size)) {
            G->capa[edge_count] = e_cost((i - 1) / grids_node + 1, sink_label);
        }
        edge_count++;
    }


    //  s->tの一連の枝から定数値を引く処理
    for (i = s2i_begin; i <= grids_node; i++) {
        G->capa[i] -= min[i];
    }
    for (i = 1; i <= grids_node; i++) {
        G->capa[i + i2t_begin - 1] -= min[i];
    }
    current_edge = depth_begin;
    for (i = 1; i <= grids_node; i++) {
        for (j = 1; j < range_size; j++) {
            G->capa[current_edge] -=min[i];
            current_edge++;
        }
    }

    free(min);
    // printf("total edge : %d\n", edge_count - 1);
    return;
}