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

double **CostMatrix(struct node_gra *net){
  struct node_gra *p, *q;
  struct node_lis *neigh;
  int nnod = CountNodes(net);
  double **cmat = allocate_d_mat(nnod, nnod);
  double links = 0, fac1=0, fac2=0, s1,s2;
  
  p = net;
  while ((p = p->next) != NULL) {
	links += (double)NodeDegree(p);
  }
  fac1 = 2/links;
  fac2 = 1/(2*links*links);

  p = net;
  while ((p = p->next) != NULL) {
	s1 = (double)NodeDegree(p);
	q = net;
	while ((q = q->next) != NULL) {
	  s2 = (double)NodeDegree(q);
	  if(cmat[q->num][p->num] != 0)
		cmat[p->num][q->num] = cmat[q->num][p->num];
	  else{
		cmat[p->num][q->num] = - fac2 * s1 * s2;
		neigh = p->neig;
		while ((neigh = neigh->next) != NULL){
		  if (neigh->node == q->num){
			cmat[p->num][q->num] += fac1;
			break;
		  }
		}
	  }
	}
  }
  return cmat;
}

double **CostMatrixWeighted(struct node_gra *net){
  struct node_gra *p, *q;
  struct node_lis *neigh;
  int nnod = CountNodes(net);
  double **cmat = allocate_d_mat(nnod, nnod);
  double links = 0, fac1=0, fac2=0, s1,s2;
  
  p = net;
  while ((p = p->next) != NULL) {
	links += (double)NodeStrength(p);
  }
  fac1 = 2/links;
  fac2 = 1/(2*links*links);

  p = net;
  while ((p = p->next) != NULL) {
	s1 = (double)NodeStrength(p);
	q = net;
	while ((q = q->next) != NULL) {
	  s2 = (double)NodeStrength(q);
	  if(cmat[q->num][p->num] != 0)
		cmat[p->num][q->num] = cmat[q->num][p->num];
	  else{
		cmat[p->num][q->num] = - fac2 * s1 * s2;
		neigh = p->neig;
		while ((neigh = neigh->next) != NULL){
		  if (neigh->node == q->num){
			cmat[p->num][q->num] += fac1*neigh->weight;
			break;
		  }
		}
	  }
	}
  }
  return cmat;
}

	
double **CostMatrixBipart(struct binet *binet){
  struct node_gra *p, *q;
  int nnod = CountNodes(binet->net1);
  double **cmat = allocate_d_mat(nnod, nnod);
  double links=0, fac1=0, fac2=0,s1=0,s2=0;
  p = binet->net2;
  
  while ((p = p->next) != NULL) {
	links = (double)NodeDegree(p);
	fac1 += links * (links - 1);
	fac2 += links;
  }
  fac1 = 1. / fac1; // 1/sum_a[m_a(m_a-1)]
  fac2 = 1. / (fac2*fac2); // 1/(sum_a[m_a])**2

  p = binet->net1;
  while ((p = p->next) != NULL) {
	s1 = (double)NodeDegree(p);
	q = binet->net1;
	while ((q = q->next) != NULL) {
	  s2 = (double)NodeDegree(q);
	  if(cmat[q->num][p->num] != 0)
		cmat[p->num][q->num] = cmat[q->num][p->num];
	  else
		cmat[p->num][q->num] = fac1 * (double)NCommonLinksBipart(p, q) - fac2 * s1 * s2;
	}
  }
  return cmat;
}


double **CostMatrixBipartWeighted(struct binet *binet){
  struct node_gra *p, *q;
  int nnod = CountNodes(binet->net1);
  struct node_lis *neigh;
  double **cmat = allocate_d_mat(nnod, nnod);
  double links, fac1, fac2, s1, s2, min;
  p = binet->net2;
  
  while ((p = p->next) != NULL) {
	links = NodeStrength(p);
	min = links;
	neigh = p->neig;
	while ((neigh = neigh->next) != NULL){
	  if (neigh->weight < min){
		min = neigh->weight;
	  }
	}
	fac1 += links * (links - min);
	fac2 += links;
  }
  fac1 = 1. / fac1; // 1/sum_a[m_a(m_a-1)]
  fac2 = 1. / (fac2*fac2); // 1/(sum_a[m_a])**2

  p = binet->net1;
  while ((p = p->next) != NULL) {
	s1 = (double)NodeStrength(p);
	q = binet->net1;
	while ((q = q->next) != NULL) {
	  s2 = (double)NodeStrength(q);
	  if(cmat[q->num][p->num] != 0)
		cmat[p->num][q->num] = cmat[q->num][p->num];
	  else
		cmat[p->num][q->num] = fac1 * (double)SumProductsOfCommonWeightsBipart(p, q) - fac2 * s1 * s2;
	}
  }
  return cmat;
}
