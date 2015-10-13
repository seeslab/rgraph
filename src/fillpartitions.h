/**
@file fillpartitions.h
@author Guilhem Doulcier
@date 2015
@license GPLv3+
@brief Convert edge list to adjacency arrays.

This file contains the functions to convert an edge list (encoded as
two arrays of unsigned int, identifying the nodes and one of double,
identifying the weights) into one AdjaArray and associated Partition
(as defined in partition.h). 

Those functions can project a bipartite edge list to perform a
bipartite SA. 

Note that those function assume that all edges are undirected, wich is
the usual case in modularity optimization.

The input arrays are usually read from a file/stream (using the io.c
file) or are coming from R arrays in the R bindings of the library.

Note: I am not still entirely satisfied with the bipartite projection
algorithm, I think you can reduce the number of times you go through
the data. Feel free to contribute.
**/

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
