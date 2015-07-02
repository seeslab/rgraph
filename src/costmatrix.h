#include "graph.h"
#include "bipartite.h"

double **CostMatrix(struct node_gra *net);
double **CostMatrixWeighted(struct node_gra *net);
double **CostMatrixBipart(struct binet *binet);
double **CostMatrixBipartWeighted(struct binet *binet);
