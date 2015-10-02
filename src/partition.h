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
