#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "bipartite.h"
#include "sannealing.h"
#include "costmatrix.h"
#include "partition.h"

/**
Normalize edges weight and node strength and store them in the
 Partition and AdjaArray objects.

This function compute the costs for Newman's modularity.

If W_ij is the adjacency matrix and W=\sum_i\sum_j W_ij the sum of its
elements. The normalized strength is:
- For an edge:  A_ij = W_ij / W
- For a node: k_i = sum_i (W_ij) / W

@param net The network to read.
@param adj The adjacency array to write.
@param part The partition to write/
**/
void
ComputeCost(struct node_gra *net,
			AdjaArray *adj, Partition *part){
  struct node_gra *node;
  struct node_lis *neig;
  double totlinks = 0;
  unsigned int i;

  node = net;
  while ((node = node->next) != NULL)
	totlinks += NodeStrength(node);

  i = 0;
  node = net;
  while((node=node->next)!=NULL){
	// Fill the AdjaArray.
	adj->idx[node->num] = i; //Start of the index of the neighbors of
							 //the focal node.
	neig = node->neig;
	while((neig=neig->next)!=NULL){
	  adj->neighbors[i] = neig->node;
	  adj->strength[i] = neig->weight/totlinks;
	  i++;
	}

	// Fill the Partition.
	part->nodes[node->num]->strength = NodeStrength(node)/totlinks;
  }
}

/**
Normalize edges weight and node strength and store them in the
 Partition and AdjaArray objects.

This function compute the costs for Guimerà's bipartite modularity.

If W is the adjacency matrix, and d_a the degree of team a.
- For an edge:  A_ij = sum_a W_{ia}W_{ja} / sum_a [d_a(d_a - epsilon_a)]
- For a node: k_i = sum_a (W_ia) / sum_a d_a

With \epsilon_a = min_{i} W_{ia}

@param net The network to read.
@param adj The adjacency array to write.
@param part The partition to write.
@param projected The projected network to read (contains the dot products sum_a[WiaWja]).
@param weighted If false, epsilon = 1.
**/
void
ComputeCostBipart(struct binet *binet,
				  AdjaArray *adj, Partition *part,
				  struct node_gra *projected,
				  unsigned int weighted){
  double epsilon;
  struct node_gra *node;
  struct node_lis *neigh;
  double links = 0;
  unsigned int i;
  double fac1,fac2;

  // Compute the normalization factors.
  node = binet->net2;
  while ((node = node->next) != NULL) {
	links = NodeStrength(node);
	if (weighted){
		epsilon = links;
		neigh = node->neig;
		while ((neigh = neigh->next) != NULL){
		  if (neigh->weight < epsilon)
			epsilon = neigh->weight;
		}
	}
	else
	  epsilon = 1.0;
	fac1 += links * (links - epsilon);
	fac2 += links;
  }
  fac1 = 1. / fac1; // 1/sum_a[m_a(m_a-epsilon)]
  fac2 = 1. / fac2; // 1/(sum_a[m_a])

  // Fill the partition with the k_i.
  node = binet->net1;
  while((node=node->next)!=NULL)
	part->nodes[node->num]->strength = NodeStrength(node)*fac2;

  // Fill the adjacency array with the A_ij.
  node = projected;
  i = 0;
  while((node=node->next)!=NULL){
	adj->idx[node->num] = i; //Start of the index of the neighbors of the focal node.
	neigh = node->neig;
	while((neigh=neigh->next)!=NULL){
	  adj->neighbors[i] = neigh->node;
	  adj->strength[i] = neigh->weight * fac1;
	  i++;
	}
  }
}
