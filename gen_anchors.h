#ifndef MATRIX_H
#define MATRIX_H
#include <string>

typedef struct matrix{
    int rows, cols;
    float **vals;
} matrix;

typedef struct {
    int *assignments;
    matrix centers;
} model;

int gen_anchors(std::string xml_dir, int num_of_clusters, int width, int height);

#endif
