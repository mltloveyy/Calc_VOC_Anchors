#define _CRT_SECURE_NO_WARNINGS
#include "gen_anchors.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include "Cascade_util.h"

#define __COMPAR_FN_T
typedef int(*__compar_fn_t)(const void*, const void*);
int anchors_data_comparator(const float **pa, const float **pb)
{
  float *a = (float *)*pa;
  float *b = (float *)*pb;
  float diff = b[0] * b[1] - a[0] * a[1];
  if (diff < 0) return 1;
  else if (diff > 0) return -1;
  return 0;
}

int count_fields(char *line)
{
  int count = 0;
  int done = 0;
  char *c;
  for (c = line; !done; ++c) {
    done = (*c == '\0');
    if (*c == ',' || done) ++count;
  }
  return count;
}

float *parse_fields(char *line, int n)
{
  float* field = (float*)calloc(n, sizeof(float));
  char *c, *p, *end;
  int count = 0;
  int done = 0;
  for (c = line, p = line; !done; ++c) {
    done = (*c == '\0');
    if (*c == ',' || done) {
      *c = '\0';
      field[count] = strtod(p, &end);
      if (p == c) field[count] = nan("");
      if (end != c && (end != c - 1 || *end != '\r')) field[count] = nan(""); //DOS file formats!
      p = c + 1;
      ++count;
    }
  }
  return field;
}

char *fgetl(FILE *fp)
{
  if (feof(fp)) return 0;
  size_t size = 512;
  char* line = (char*)malloc(size * sizeof(char));
  if (!fgets(line, size, fp)) {
    free(line);
    return 0;
  }

  size_t curr = strlen(line);

  while ((line[curr - 1] != '\n') && !feof(fp)) {
    if (curr == size - 1) {
      size *= 2;
      line = (char*)realloc(line, size * sizeof(char));
    }
    size_t readsize = size - curr;
    if (readsize > INT_MAX) readsize = INT_MAX - 1;
    fgets(&line[curr], readsize, fp);
    curr = strlen(line);
  }
  if (curr >= 2)
    if (line[curr - 2] == 0x0d) line[curr - 2] = 0x00;

  if (curr >= 1)
    if (line[curr - 1] == 0x0a) line[curr - 1] = 0x00;

  return line;
}

void top_k(float *a, int n, int k, int *index)
{
  int i, j;
  for (j = 0; j < k; ++j) index[j] = -1;
  for (i = 0; i < n; ++i) {
    int curr = i;
    for (j = 0; j < k; ++j) {
      if ((index[j] < 0) || a[curr] > a[index[j]]) {
        int swap = curr;
        curr = index[j];
        index[j] = swap;
      }
    }
  }
}

void file_error(char *s)
{
  fprintf(stderr, "Couldn't open file: %s\n", s);
  exit(EXIT_FAILURE);
}

void free_matrix(matrix m)
{
    int i;
    for(i = 0; i < m.rows; ++i) free(m.vals[i]);
    free(m.vals);
}

float matrix_topk_accuracy(matrix truth, matrix guess, int k)
{
    int* indexes = (int*)calloc(k, sizeof(int));
    int n = truth.cols;
    int i,j;
    int correct = 0;
    for(i = 0; i < truth.rows; ++i){
        top_k(guess.vals[i], n, k, indexes);
        for(j = 0; j < k; ++j){
            int class_id = indexes[j];
            if(truth.vals[i][class_id]){
                ++correct;
                break;
            }
        }
    }
    free(indexes);
    return (float)correct/truth.rows;
}

void scale_matrix(matrix m, float scale)
{
    int i,j;
    for(i = 0; i < m.rows; ++i){
        for(j = 0; j < m.cols; ++j){
            m.vals[i][j] *= scale;
        }
    }
}

matrix resize_matrix(matrix m, int size)
{
    int i;
    if (m.rows == size) return m;
    if (m.rows < size) {
        m.vals = (float**)realloc(m.vals, size * sizeof(float*));
        for (i = m.rows; i < size; ++i) {
            m.vals[i] = (float*)calloc(m.cols, sizeof(float));
        }
    } else if (m.rows > size) {
        for (i = size; i < m.rows; ++i) {
            free(m.vals[i]);
        }
        m.vals = (float**)realloc(m.vals, size * sizeof(float*));
    }
    m.rows = size;
    return m;
}

void matrix_add_matrix(matrix from, matrix to)
{
    assert(from.rows == to.rows && from.cols == to.cols);
    int i,j;
    for(i = 0; i < from.rows; ++i){
        for(j = 0; j < from.cols; ++j){
            to.vals[i][j] += from.vals[i][j];
        }
    }
}

matrix make_matrix(int rows, int cols)
{
    int i;
    matrix m;
    m.rows = rows;
    m.cols = cols;
    m.vals = (float**)calloc(m.rows, sizeof(float*));
    for(i = 0; i < m.rows; ++i){
        m.vals[i] = (float*)calloc(m.cols, sizeof(float));
    }
    return m;
}

matrix hold_out_matrix(matrix *m, int n)
{
    int i;
    matrix h;
    h.rows = n;
    h.cols = m->cols;
    h.vals = (float**)calloc(h.rows, sizeof(float*));
    for(i = 0; i < n; ++i){
        int index = rand()%m->rows;
        h.vals[i] = m->vals[index];
        m->vals[index] = m->vals[--(m->rows)];
    }
    return h;
}

float *pop_column(matrix *m, int c)
{
    float* col = (float*)calloc(m->rows, sizeof(float));
    int i, j;
    for(i = 0; i < m->rows; ++i){
        col[i] = m->vals[i][c];
        for(j = c; j < m->cols-1; ++j){
            m->vals[i][j] = m->vals[i][j+1];
        }
    }
    --m->cols;
    return col;
}

matrix csv_to_matrix(char *filename)
{
    FILE *fp = fopen(filename, "r");
    if(!fp) file_error(filename);

    matrix m;
    m.cols = -1;

    char *line;

    int n = 0;
    int size = 1024;
    m.vals = (float**)calloc(size, sizeof(float*));
    while((line = fgetl(fp))){
        if(m.cols == -1) m.cols = count_fields(line);
        if(n == size){
            size *= 2;
            m.vals = (float**)realloc(m.vals, size * sizeof(float*));
        }
        m.vals[n] = parse_fields(line, m.cols);
        free(line);
        ++n;
    }
    m.vals = (float**)realloc(m.vals, n * sizeof(float*));
    m.rows = n;
    return m;
}

void matrix_to_csv(matrix m)
{
    int i, j;

    for(i = 0; i < m.rows; ++i){
        for(j = 0; j < m.cols; ++j){
            if(j > 0) printf(",");
            printf("%.17g", m.vals[i][j]);
        }
        printf("\n");
    }
}

void print_matrix(matrix m)
{
    int i, j;
    printf("%d X %d Matrix:\n",m.rows, m.cols);
    printf(" __");
    for(j = 0; j < 16*m.cols-1; ++j) printf(" ");
    printf("__ \n");

    printf("|  ");
    for(j = 0; j < 16*m.cols-1; ++j) printf(" ");
    printf("  |\n");

    for(i = 0; i < m.rows; ++i){
        printf("|  ");
        for(j = 0; j < m.cols; ++j){
            printf("%15.7f ", m.vals[i][j]);
        }
        printf(" |\n");
    }
    printf("|__");
    for(j = 0; j < 16*m.cols-1; ++j) printf(" ");
    printf("__|\n");
}


matrix make_matrix(int rows, int cols);

void copy(float *x, float *y, int n);
float dist(float *x, float *y, int n);
int *sample(int n);

int closest_center(float *datum, matrix centers)
{
    int j;
    int best = 0;
    float best_dist = dist(datum, centers.vals[best], centers.cols);
    for (j = 0; j < centers.rows; ++j) {
        float new_dist = dist(datum, centers.vals[j], centers.cols);
        if (new_dist < best_dist) {
            best_dist = new_dist;
            best = j;
        }
    }
    return best;
}

float dist_to_closest_center(float *datum, matrix centers)
{
    int ci = closest_center(datum, centers);
    return dist(datum, centers.vals[ci], centers.cols);
}

int kmeans_expectation(matrix data, int *assignments, matrix centers)
{
    int i;
    int converged = 1;
    for (i = 0; i < data.rows; ++i) {
        int closest = closest_center(data.vals[i], centers);
        if (closest != assignments[i]) converged = 0;
        assignments[i] = closest;
    }
    return converged;
}

void kmeans_maximization(matrix data, int *assignments, matrix centers)
{
    matrix old_centers = make_matrix(centers.rows, centers.cols);

    int i, j;
    int *counts = (int*)calloc(centers.rows, sizeof(int));
    for (i = 0; i < centers.rows; ++i) {
        for (j = 0; j < centers.cols; ++j) {
            old_centers.vals[i][j] = centers.vals[i][j];
            centers.vals[i][j] = 0;
        }
    }
    for (i = 0; i < data.rows; ++i) {
        ++counts[assignments[i]];
        for (j = 0; j < data.cols; ++j) {
            centers.vals[assignments[i]][j] += data.vals[i][j];
        }
    }
    for (i = 0; i < centers.rows; ++i) {
        if (counts[i]) {
            for (j = 0; j < centers.cols; ++j) {
                centers.vals[i][j] /= counts[i];
            }
        }
    }

    for (i = 0; i < centers.rows; ++i) {
        for (j = 0; j < centers.cols; ++j) {
            if(centers.vals[i][j] == 0) centers.vals[i][j] = old_centers.vals[i][j];
        }
    }
    free(counts);
    free_matrix(old_centers);
}



void random_centers(matrix data, matrix centers) {
    int i;
    int *s = sample(data.rows);
    for (i = 0; i < centers.rows; ++i) {
        copy(data.vals[s[i]], centers.vals[i], data.cols);
    }
    free(s);
}

int *sample(int n)
{
    int i;
    int* s = (int*)calloc(n, sizeof(int));
    for (i = 0; i < n; ++i) s[i] = i;
    for (i = n - 1; i >= 0; --i) {
        int swap = s[i];
        int index = rand() % (i + 1);
        s[i] = s[index];
        s[index] = swap;
    }
    return s;
}

float dist(float *x, float *y, int n)
{
    //printf(" x0 = %f, x1 = %f, y0 = %f, y1 = %f \n", x[0], x[1], y[0], y[1]);
    float mw = (x[0] < y[0]) ? x[0] : y[0];
    float mh = (x[1] < y[1]) ? x[1] : y[1];
    float inter = mw*mh;
    float sum = x[0] * x[1] + y[0] * y[1];
    float un = sum - inter;
    float iou = inter / un;
    return 1 - iou;
}

void copy(float *x, float *y, int n)
{
    int i;
    for (i = 0; i < n; ++i) y[i] = x[i];
}

model do_kmeans(matrix data, int k)
{
    matrix centers = make_matrix(k, data.cols);
    int* assignments = (int*)calloc(data.rows, sizeof(int));
    //smart_centers(data, centers);
    random_centers(data, centers);  // IoU = 67.31% after kmeans

    /*
    // IoU = 63.29%, anchors = 10,13,  16,30,  33,23,  30,61,  62,45,  59,119,  116,90,  156,198,  373,326
    centers.vals[0][0] = 10; centers.vals[0][1] = 13;
    centers.vals[1][0] = 16; centers.vals[1][1] = 30;
    centers.vals[2][0] = 33; centers.vals[2][1] = 23;
    centers.vals[3][0] = 30; centers.vals[3][1] = 61;
    centers.vals[4][0] = 62; centers.vals[4][1] = 45;
    centers.vals[5][0] = 59; centers.vals[5][1] = 119;
    centers.vals[6][0] = 116; centers.vals[6][1] = 90;
    centers.vals[7][0] = 156; centers.vals[7][1] = 198;
    centers.vals[8][0] = 373; centers.vals[8][1] = 326;
    */

    // range centers [min - max] using exp graph or Pyth example
    if (k == 1) kmeans_maximization(data, assignments, centers);
    int i;
    for(i = 0; i < 1000 && !kmeans_expectation(data, assignments, centers); ++i) {
        kmeans_maximization(data, assignments, centers);
    }
    printf("\n iterations = %d \n", i);
    model m;
    m.assignments = assignments;
    m.centers = centers;
    return m;
}

int gen_anchors(std::string xml_dir, int num_of_clusters, int width, int height)
{
  int files_number;
  std::vector<std::string> files_name;
  readFiles(xml_dir + "*.xml", files_name, &files_number);
  int number_of_boxes = 0;
  float* rel_width_height_array = (float*)calloc(1000, sizeof(float));

  for (int i = 0; i < files_name.size(); i++)
  {
    std::string file_name = xml_dir + files_name[i];
    xmlReadWrite img_info;
    readXmlFile(file_name, img_info, "");
    int ori_width = img_info.img_width;
    int ori_height = img_info.img_hight;
    for (int j = 0; j < img_info.xml_object_vec.size(); j++)
    {
      cv::Rect box_i = img_info.xml_object_vec[j].gt_boxes;

      number_of_boxes++;
      rel_width_height_array = (float*)realloc(rel_width_height_array, 2 * number_of_boxes * sizeof(float));

      rel_width_height_array[number_of_boxes * 2 - 2] = 1.0 * box_i.width / ori_width * width;
      rel_width_height_array[number_of_boxes * 2 - 1] = 1.0 * box_i.height / ori_height * height;
      std::printf("\r loaded  image: %d  box: %d", i + 1, number_of_boxes);
    }
  }

  std::printf("\n all loaded. \n");
  std::printf("\n calculating k-means++ ...");

  matrix boxes_data;
  model anchors_data;
  boxes_data = make_matrix(number_of_boxes, 2);

  std::printf("\n");
  for (int i = 0; i < number_of_boxes; ++i) {
    boxes_data.vals[i][0] = rel_width_height_array[i * 2];
    boxes_data.vals[i][1] = rel_width_height_array[i * 2 + 1];
  }

  // K-means
  anchors_data = do_kmeans(boxes_data, num_of_clusters);

  std::qsort((void*)anchors_data.centers.vals, num_of_clusters, 2 * sizeof(float), (__compar_fn_t)anchors_data_comparator);

  std::printf("\n");
  float avg_iou = 0;
  for (int i = 0; i < number_of_boxes; ++i) {
    float box_w = rel_width_height_array[i * 2]; //points->data.fl[i * 2];
    float box_h = rel_width_height_array[i * 2 + 1]; //points->data.fl[i * 2 + 1];

    int cluster_idx = 0;
    float min_dist = FLT_MAX;
    float best_iou = 0;
    for (int j = 0; j < num_of_clusters; ++j) {
      float anchor_w = anchors_data.centers.vals[j][0];   // centers->data.fl[j * 2];
      float anchor_h = anchors_data.centers.vals[j][1];   // centers->data.fl[j * 2 + 1];
      float min_w = (box_w < anchor_w) ? box_w : anchor_w;
      float min_h = (box_h < anchor_h) ? box_h : anchor_h;
      float box_intersect = min_w*min_h;
      float box_union = box_w*box_h + anchor_w*anchor_h - box_intersect;
      float iou = box_intersect / box_union;
      float distance = 1 - iou;
      if (distance < min_dist) {
        min_dist = distance;
        cluster_idx = j;
        best_iou = iou;
      }
    }

    float anchor_w = anchors_data.centers.vals[cluster_idx][0]; //centers->data.fl[cluster_idx * 2];
    float anchor_h = anchors_data.centers.vals[cluster_idx][1]; //centers->data.fl[cluster_idx * 2 + 1];
    if (best_iou > 1 || best_iou < 0) { // || box_w > width || box_h > height) {
      std::printf(" Wrong label: i = %d, box_w = %f, box_h = %f, anchor_w = %f, anchor_h = %f, iou = %f \n",
        i, box_w, box_h, anchor_w, anchor_h, best_iou);
    }
    else avg_iou += best_iou;
  }

  std::printf(" Anchor = \n");
  for (int i = 0; i < num_of_clusters; ++i) {
    float anchor_w = anchors_data.centers.vals[i][0]; //centers->data.fl[i * 2];
    float anchor_h = anchors_data.centers.vals[i][1]; //centers->data.fl[i * 2 + 1];
    std::printf("%4.0f, %4.0f, %4.0f, %4.0f\n", -anchor_w / 2, -anchor_h / 2, anchor_w / 2, anchor_h / 2);
  }

  avg_iou = 100 * avg_iou / number_of_boxes;
  std::printf("\n avg IoU = %2.2f %% \n", avg_iou);
  std::free(rel_width_height_array);
}
