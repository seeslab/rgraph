/**
@file partitions.h
@author Guilhem Doulcier
@date 2015
@license GPLv3+
@brief Data structure optimized for modularity optimization by SA.

This file contains the definition and useful methods for the
Partition and AdjaArray (Adjacency Array) data structure which are at
the core of this implementation of the SA algorithm.

We need an efficient way to:
1. Compute the modularity variation due to changing the module of a vertex. 
2. Actually change the module of a vertex.
3. Compute the modularity variation due to the merging of two modules. 
4. Actually merge two modules.

The adjacency array and the partition as a list of doubly linked lists
is an efficient solution to this problem, even if the graph is
relatively static.

To loop through:
- All nodes: use the Partition->nodes array.
- Neighbors of node i: use the AdjaArray->neighbors array between id[i] and id[i+1]-1.
- Nodes member of a module j: use the doubly linked list partition->module[j]. 
**/

#ifndef PARTITION_H__
#define PARTITION_H__
#include <gsl/gsl_rng.h>

/**
Adjacency Array.

This structure is made to store a static undirected graph with N nodes
and E edges.

Given a node i, its neighbors are the nodes indexed by
neighbors[idx[i]] through neighbors[idx[i]+1].

- idx is an array of size N.
- neighbors and strength are arrays of size 2E+1 and 2E.
**/
typedef struct AdjaArray {
  unsigned int N; //! Number of nodes.
  unsigned int E; //! Number of edges.
  unsigned int *idx; //!< idx[i] is the index of the first neighbor of i in neighbors.
  unsigned int *neighbors; //!< neighbors[idx[i]] to neighbors[idx[i]+1] are the neighbors of i.
  double *strength; //!< strength of the corresponding edge in neighbors.
} AdjaArray;

AdjaArray * CreateAdjaArray(unsigned int N, unsigned int E);
void FreeAdjaArray(AdjaArray *adj);

/**
Node.

A doubly linked list of nodes.
   **/
typedef struct Node{
  unsigned int id; //! Node id (to use in the adjacency array for instance).
  double strength; //! (normalized) strength.
  unsigned int module; //! module id.
  struct Node *prev; //!< the previous node in the module.
  struct Node *next; //!< the next node in the module.
} Node;

/**
Module.

This structure is essentially a pointer toward the doubly link list of
the nodes inside the module, with a few variables for book-keeping.
**/
typedef struct Module {
  unsigned int id;
  unsigned int size; //!< Number of node currently in the module.
  double strength; //!< Sum of the strength of nodes in the module.
  Node *first; //!< First node of the module.
  Node *last; //!< Last node of the module.
} Module;

/**
Partition
**/
typedef struct Partition {
  unsigned int N; //!< Number of nodes.
  unsigned int M; //!< Number of modules.
  unsigned int nempty; //!< Number of empty modules.
  Module **modules; //! An array of modules.
  Node **nodes; //! An array of pointers to nodes.
} Partition;


Partition * CreatePartition(unsigned int N, unsigned int M);
Partition * CopyPartitionStruct(Partition *part);
void FreePartition(Partition *part);

void AssignNodesToModules(Partition *part, 	gsl_rng *gen);
double PartitionModularity(Partition *part, AdjaArray *adj, int diagonal_term);
void PartitionRolesMetrics(Partition *part, AdjaArray *adj,  double connectivity[], double participartion[]);
void CompressPartition(Partition *part);
int GetRole(double P, double z);
#endif
