#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bmp.h"
#include "rangeswap.h"
#include "graph.h"
#include "ford_fulkerson.h"

#define _OUTPUT_T_ 1     // BK-maxflow後のtの状態をファイルに出力 0:出力しない 1:出力する
#define _OUTPUT_INFO_ 0     // デバッグ情報出力 0:出力しない 1:出力する
#define _OUTPUT_GRAPH_ 0    // グラフ情報出力  0:出力しない 1:出力する
#define _OUTPUT_PROGRESS_ 1 // 処理過程ファイル出力 0:出力しない 1:出力する
#define _RUN_FIRST_ONLY_ 0 //1度目の移動で終了(デバッグ用)

int main(int argc, char *argv[]) {
    int i, j, k, node, edge, grids_node, flag, alpha, beta, swap_node_size;
    int scale, label_max, grids_edge, count, last_move, num_of_moves;
    int *I, *t, *label, *newlabel, *label_index;
    // I->入力画像の輝度, t->2値変数, label->ラベル付け
    int label_size = 8;
    int range_size = 4;
    double decreace, prev_energy, T = 255;
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
#if _OUTPUT_PROGRESS_
    int l = 0;
    char pf[100];
    system("rm output/*.bmp &> /dev/null");
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
    if (range_size < 2) {
        fprintf(stderr, "Error! Range size == %d \n", range_size);
        exit (EXIT_FAILURE);
    }

    printf("----------------------------------------------\n");
    printf("input_file: %s\n", argv[1]);
    printf("output_file: %s\n", output_file);
    printf("label_size: %d\n", label_size);
    if(h(2) > 2) printf("Vpq(fp, fq) = (fp - fq)^2\n");
    else printf("Vpq(fp, fq) = |fp - fq|\n");

    ReadBmp(argv[1], &image);
    ReadBmp(argv[1], &output);

    grids_node = image.height * image.width;


    if ((I = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->I]\n");
        exit(EXIT_FAILURE);
    }

    printf("height %ld, width %ld\n", image.height, image.width);
    for (i = 0; i <  image.height; i++) {
        for (j = 0; j < image.width; j++) {
            if(image.data[i][j].r / scale > label_max) I[i * image.width + j + 1] = label_max;
            else I[i * image.width + j + 1] = image.data[i][j].r / scale;
        }
    }

    // エネルギー計算用一層グラフ作成
    node = grids_node + 2;
    edge = (image.height - 1) * image.width + image.height * (image.width - 1) + 2 * grids_node;
    newGraph(&Ge, node, edge);

    // set_single_edge(&Ge, image.height, image.width);
    set_single_edges(&Ge, image.height, image.width);
    initAdjList(&Ge);

    if ((label = (int *) malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [label = malloc()]\n");
        return (EXIT_FAILURE);
    }
    if ((newlabel = (int *) malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [label = malloc()]\n");
        return (EXIT_FAILURE);
    }
    if ((label_index = (int *) malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [label_index = malloc()]\n");
        return (EXIT_FAILURE);
    }
    
    
    // 輝度から初期ラベル設定
    for (i = 1; i <= grids_node ; i++) label[i] = I[i];
    cpyarray(newlabel, label, grids_node);
    prev_energy = energy(&Ge, label, I, T);
    printf("Energy (before): %lf\n", prev_energy);


#if _OUTPUT_T_
    fprintf(fp, "Energy (before): %lf\n", prev_energy);
    fprintf(fp, "init_label:\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", label[i]);
        if(i % image.width == 0) fprintf(fp, "\n");
        if(i % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
    }

#endif

    num_of_moves = label_size - range_size + 1;
    last_move = num_of_moves;
    decreace = 0;
    flag = 0;
    start = clock();
    do {
        prev_energy = energy(&Ge, label, I, T);
        for(i = 0; i < num_of_moves; i++) {
            alpha = i;
            beta = alpha + range_size - 1;
            printf("alpha: %d beta: %d\n", alpha, beta);
            if (last_move == i) {
                flag = 1;
                break;
            }
            swap_node_size = make_label_index(&Ge, label, label_index, alpha, beta);
           
            // for (j = 1; j <= Ge.n - 2; j++) {
            //     // printf("t[%d] : %d\n", i, t[i]);
            //     printf("%d ", isin_array(label_index, j, swap_node_size) ? 1 : 0);
            //     if(j % image.width == 0) printf("\n");
            //     if(j % (grids_node) == 0) printf("-------------------------------------\n");
            // }

            if (swap_node_size  < 2) continue;

#if _OUTPUT_T_
            fprintf(fp, "\n-------------------------------------\n");
            fprintf(fp, "alpha: %d beta: %d\n", alpha, beta);

            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", isin_array(label_index, j, swap_node_size) ? 1 : 0);
                if(j % image.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
            fprintf(fp, "label: \n");
            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", label[j]);
                if(j % image.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }      
#endif
        
            
            node = image.height * image.width * (range_size - 1) + 2;
            grids_edge = (image.height - 1) * image.width + image.height * (image.width - 1);
            edge =  2 * ((range_size - 1) * ((range_size - 1) * grids_edge + 2 * grids_node) + grids_node * ((range_size - 1) - 1));
            
            newGraph(&G, node, edge);
            set_edge(&G, image.height, image.width, alpha, beta, label_size, label, label_index, swap_node_size, I);
            initAdjList(&G);
            
            if ((f = (double *) malloc(sizeof(double) * (G.m + 1))) == NULL) {
                fprintf(stderr, "main(): ERROR [f = malloc()]\n");
                return (EXIT_FAILURE);
            }
            if ((t = (int *) malloc(sizeof(int) * (G.n + 1))) == NULL) {
                fprintf(stderr, "main(): ERROR [t = malloc()]\n");
                return (EXIT_FAILURE);
            }

            for (j = 0; j < G.m + 1 ; j++) f[j] = 0;
            for (j = 0; j < G.n + 1 ; j++) t[j] = 0;

            boykov_kolmogorov(G, f, t);
            for (j = 1; j <= Ge.n - 2; j++) {
                if (isin_array(label_index, j, Ge.n - 2)) {
                    k = j;
                    count = 0;
                    while (t[k] == 1 && k <= (range_size - 1) * grids_node) {
                        k += grids_node;
                        count++;
                    }
                    newlabel[j] = count + alpha;
                } else newlabel[j] = label[j];
            }

            if (cmparray(newlabel, label, grids_node)) {
                last_move = alpha;
                cpyarray(label, newlabel, grids_node);
                printf("Energy : %lf\n", energy(&Ge, label, I, T));
            }

#if _OUTPUT_T_
            fprintf(fp, "t: \n");
            for (j = 1; j <= G.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", t[j]);
                if(j % image.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
            fprintf(fp, "label: \n");
            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", newlabel[j]);
                if(j % image.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }

            
#endif
#if _OUTPUT_PROGRESS_
            for (j = 0; j <  image.height; j++) {
                for (k = 0; k < image.width; k++) {
                    output.data[j][k].r = label[j * image.width + k + 1] * scale;
                    output.data[j][k].g = output.data[j][k].r;
                    output.data[j][k].b = output.data[j][k].r;
                }
            }
            sprintf(pf, "output/image_%04d.bmp", l);
            WriteBmp(pf, &output);
            l++;
#endif
            // showGraph(&G);
            free(f);
            free(t);
            delGraph(&G);

#if _RUN_FIRST_ONLY_
            flag = 1;
            break;
 #endif
        }
        if (flag) break;
        decreace = prev_energy - energy(&Ge, label, I, T);
        printf("Energy : %lf\n", energy(&Ge, label, I, T));
    } while (decreace > 0);
    
#if _OUTPUT_T_
    fprintf(fp, "result:\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        fprintf(fp, "%d ", label[i]);
        if(i % image.width == 0) fprintf(fp, "\n");
    }
#endif

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

    
    delGraph(&Ge);
    free(I);
    free(label);
    free(newlabel);
    free(label_index);
    printf("----------------------------------------------\n");
    return 0;
    
}
