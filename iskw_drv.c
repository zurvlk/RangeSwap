#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bmp.h"
#include "ishikawa.h"
#include "graph.h"
#include "ford_fulkerson.h"

#define _OUTPUT_T_ 0     // BK-maxflow後のtの状態をファイルに出力 0:出力しない 1:出力する
#define _OUTPUT_INFO_ 0     // デバッグ情報出力 0:出力しない 1:出力する
#define _OUTPUT_GRAPH_ 0    // グラフ情報出力 　0:出力しない 1:出力する

int main(int argc, char *argv[]) {
    int i, j, k, node, edge, grids_node;
    int scale, label_max, grids_edge, count;
    int *I, *t, *label;
    // I->入力画像の輝度, t->2値変数, label->ラベル付け
    int label_size = 16;
    double  T = 255;
    double *f;
    char output_file[100];
    clock_t start;
    img image, output;
    // Ge:エネルギー計算用
    Graph G, Ge;
    

#if _OUTPUT_INFO_
    double maxflow;
#endif
#if _OUTPUT_T_
    FILE *fp;
    fp = fopen("log/t.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "cannot open file[t.txt]\n");
        exit (EXIT_FAILURE);
    }
#endif

    if (argc != 2 && argc != 4) {
        printf("Usage: %s <input_file> <output_file(option)> <label_size(option)>\n", argv[0]);
        return 1;
    }

    if (argc == 2) strcpy(output_file, "/dev/null");
    else strcpy(output_file, argv[2]);
    if (argc == 4) label_size = atoi(argv[3]);

    label_max = label_size - 1;
    scale = 256 / label_size;
    printf("----------------------------------------------\n");
    printf("input_file: %s\n", argv[1]);
    printf("output_file: %s\n", output_file);
    printf("label_size: %d\n", label_size);
    if(h(2) > 2) printf("Vpq(fp, fq) = (fp - fq)^2\n");
    else printf("Vpq(fp, fq) = |fp - fq|\n");

    ReadBmp(argv[1], &image);
    ReadBmp(argv[1], &output);

    grids_node = image.height * image.width;


    if ((I = (int *)malloc(sizeof(int) * (image.height * image.width + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->I]\n");
        exit(EXIT_FAILURE);
    }

    printf("height %ld, width %ld\n", image.height, image.width);
    for (i = 0; i <  image.height; i++) {
        for (j = 0; j < image.width; j++) {
            I[i * image.width + j + 1] = image.data[i][j].r / scale;
        }
    }

    // エネルギー計算用一層グラフ作成
    node = image.height * image.width + 2;
    edge = (image.height - 1) * image.width + image.height * (image.width - 1) + 2 * image.height * image.width;
    newGraph(&Ge, node, edge);

    
    // set_single_edge(&Ge, image.height, image.width);
    set_single_edges(&Ge, image.height, image.width);
    initAdjList(&Ge);

    // Ishikawaのアルゴリズム用グラフのnode,edge数計算
    node = image.height * image.width * label_max + 2;
    grids_edge = (image.height - 1) * image.width + image.height * (image.width - 1);
    
    // edge =  2 * (label_max * label_max * grids_edge + grids_node * (label_max - 1) + 2 * grids_node * label_max);
    edge =  2 * (label_max * (label_max * grids_edge + 2 * grids_node) + grids_node * (label_max - 1));

    // グラフ初期設定
    
    // printf("node: %d edge : %d\n", node, edge);
    newGraph(&G, node, edge);
    set_all_edge(&G, image.height, image.width, label_max, I);

#if _OUTPUT_GRAPH_
    printf("---Graph information---\n");
    showGraph(&G);
#endif

    initAdjList(&G);
    // printf("complete initAdjList()\n");

    if ((f = (double *) malloc(sizeof(double) * (G.m + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [f = malloc()]\n");
        return (EXIT_FAILURE);
    }
    if ((t = (int *) malloc(sizeof(int) * (G.n + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [t = malloc()]\n");
        return (EXIT_FAILURE);
    }
    if ((label = (int *) malloc(sizeof(int) * (image.height * image.width + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [label = malloc()]\n");
        return (EXIT_FAILURE);
    }

    // 輝度から初期ラベル設定
    for (i = 1; i <= image.width * image.height ; i++) label[i] = 0;
    printf("Energy (before): %lf\n", energy(&Ge, label, I, T));


#if _OUTPUT_T_
    fprintf(fp, "Energy (before): %lf\n", energy(&Ge, label, I, T));
    fprintf(fp, "init_label:\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", label[i]);
        if(i % image.width == 0) fprintf(fp, "\n");
        if(i % (image.width * image.height) == 0) fprintf(fp, "-------------------------------------\n");
    }
    int sum = 0;
    for (i = 1; i <= Ge.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", I[i]);
        sum += I[i];
        if(i % image.width == 0) fprintf(fp, "\n");
        if(i % (image.width * image.height) == 0) fprintf(fp, "-------------------------------------\n");
    }
    fprintf(fp, "Data %d\n", sum);
#endif

    for (i = 0; i < G.m + 1 ; i++) f[i] = 0;
    for (i = 0; i < G.n + 1 ; i++) t[i] = 0;

    start = clock();
    /*
    printf("max flow : %lf\n", boykov_kolmogorov(G, f, t));
    /*/
    boykov_kolmogorov(G, f, t);
    // */

#if _OUTPUT_T_
    for (i = 1; i <= G.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", t[i]);
        if(i % image.width == 0) fprintf(fp, "\n");
        if(i % (image.width * image.height) == 0) fprintf(fp, "-------------------------------------\n");
    }
    fprintf(fp, "result:\n");
#endif

    for (i = 1; i <= Ge.n - 2; i++) {
        k = i;
        count = 0;
        while (t[k] == 1 && k <= label_max * image.width * image.height) {
            k += image.height * image.width;
            count++;
        }
        label[i] = count;
        
#if _OUTPUT_T_
        fprintf(fp, "%d ", label[i]);
        if(i % image.width == 0) fprintf(fp, "\n");
#endif
    }

    printf("Energy (after): %lf\n", energy(&Ge, label, I, T));
    printf("Run time[%.2lf]\n", (double) (clock() - start) / CLOCKS_PER_SEC);

    // output to bitmap file
    for (i = 0; i <  image.height; i++) {
        for (j = 0; j < image.width; j++) {
            output.data[i][j].r = label[i * image.width + j + 1] * scale;
            output.data[i][j].g = output.data[i][j].r;
            output.data[i][j].b = output.data[i][j].r;
        }
    }
    WriteBmp(output_file, &output);

#if _OUTPUT_T_
    fprintf(fp, "Energy (after): %lf\n", energy(&Ge, label, I, T));
    fclose(fp);
#endif

    delGraph(&G);
    delGraph(&Ge);
    free(I);
    free(f);
    free(t);
    free(label);
    printf("----------------------------------------------\n");
    return 0;
    
}