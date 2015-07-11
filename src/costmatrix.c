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

 **/
void
ComputeCost(struct node_gra *net,
			AdjaArray *adj, Partition *part){
  struct node_gra *node;
  struct node_lis *neig;
  double links = 0;
  unsigned int i;
  double fac;


  node = net;
  while ((node = node->next) != NULL) 
	links += NodeStrength(node);

  i = 0;
  node = net;
  while((node=node->next)!=NULL){
	adj->idx[node->num] = i; //Start of the index of the neighbors of the focal node.
	neig = node->neig;
	while((neig=neig->next)!=NULL){
	  adj->neighbors[i] = neig->node;
	  adj->strength[i] = neig->weight/links;
	  i++;
	}	
	part->nodes[node->num]->strength = NodeStrength(node)/links;
  }
}

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

  node = binet->net1;
  while((node=node->next)!=NULL)
	part->nodes[node->num]->strength = NodeStrength(node)*fac2;

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

