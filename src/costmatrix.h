#include "graph.h"
#include "bipartite.h"
#include "partition.h"

void ComputeCost(struct node_gra *net, AdjaArray *adj, Partition *part);
void ComputeCostBipart(struct binet *binet, AdjaArray *adj, Partition *part, struct node_gra *projected, unsigned int weighted);
