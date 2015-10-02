#ifndef FILLPART_H__
#define FILLPART_H__
#include "partition.h"

typedef struct Edge {
  unsigned int node1;
  unsigned int node2;
  double strength;
} Edge;

static int
EdgeCompare(const void *p1, const void *p2);

int
EdgeListToAdjaArray(int *nd_in, int *nd_out, double *weight,
					          AdjaArray *adj, Partition *part, int normalize);
unsigned int
ProjectBipartEdgeList(unsigned int *nd_in, unsigned int *nd_out, double *weights, int E,
                      Partition **part_p, AdjaArray **adj_p );


void
AssignNodesToModules(Partition *part, gsl_rng *gen);
#endif
