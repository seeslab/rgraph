/*
  models.c
  $LastChangedDate$
  $Revision$
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "models.h"

/*
  ---------------------------------------------------------------------
  Create an empty graph
  ---------------------------------------------------------------------
*/
struct node_gra *
EmptyGraph(int S)
{
  int node1;
  struct node_gra *root = NULL, *last = NULL;
  char label[MAX_LABEL_LENGTH];

  /* Create header */
  last = root = CreateHeaderGraph();
  
  /* Create nodes */
  for (node1=0; node1<S; node1++) {
    sprintf(&label[0], "%d", node1+1);
    last = CreateNodeGraph(last, &label[0]);
  }

  /* Done */
  return root;
}

/*
  ---------------------------------------------------------------------
  Create an Erdos-Renyi random graph
  ---------------------------------------------------------------------
*/
struct node_gra *
ERGraph(int S, double p, gsl_rng *gen)
{
  int node1, node2;
  struct node_gra **nodeList = NULL;
  struct node_gra *root = NULL, *last = NULL;
  char label[MAX_LABEL_LENGTH];

  /* Create header */
  last = root = CreateHeaderGraph();
  
  /* Create nodes */
  nodeList = (struct node_gra **)calloc(S, sizeof(struct node_gra *));
  for (node1=0; node1<S; node1++) {
    sprintf(&label[0], "%d", node1+1);
    last = nodeList[node1] = CreateNodeGraph(last, &label[0]);
  }

  /* Create the links */
  for (node1=0; node1<S; node1++) {
    for (node2=node1+1; node2<S; node2++) {
      if (gsl_rng_uniform(gen) < p) {
	AddAdjacency(nodeList[node1], nodeList[node2], 0, 0, 0, 0);
	AddAdjacency(nodeList[node2], nodeList[node1], 0, 0, 0, 0);
      }
    }
  }

  /* Done */
  free(nodeList);
  return root;
}

/*
  ---------------------------------------------------------------------
  Create a preferential attachment (aka Barabasi-Albert) network
  ---------------------------------------------------------------------
*/
struct node_gra *
PAGraph(int S, int m, gsl_rng *gen)
{
  int node1, node2, n;
  struct node_gra **nodeList=NULL;
  struct node_gra *root=NULL, *last=NULL;
  char label[MAX_LABEL_LENGTH];
  int *paList=NULL, norm=0;

  /* The list for selecting nodes preferentially */
  paList = (int *)calloc(2 * S * m, sizeof(int));

  /* Initial m fully connected nodes */
  nodeList = (struct node_gra **)calloc(S, sizeof(struct node_gra *));
  last = root = CreateHeaderGraph();
  for (node1=0; node1<m; node1++) {
    sprintf(&label[0], "%d", node1+1);
    last = nodeList[node1] = CreateNodeGraph(last, &label[0]);
    for (node2=0; node2<node1; node2++) {
      AddAdjacency(nodeList[node1], nodeList[node2], 0, 0, 0, 0);
      AddAdjacency(nodeList[node2], nodeList[node1], 0, 0, 0, 0);
      paList[norm++]=node1;
      paList[norm++]=node2;
    }
  }
  
  /* Create remaining nodes and links */
  for (node1=m; node1<S; node1++) {
    sprintf(&label[0], "%d", node1+1);
    last = nodeList[node1] = CreateNodeGraph(last, &label[0]);
    for (n=0; n<m; n++) {
      do {
	/* preferential attachment */
	node2 = paList[(int)(gsl_rng_uniform(gen) * norm)];
      }
      while (IsThereLink(nodeList[node1], nodeList[node2]) == 1);
      AddAdjacency(nodeList[node1], nodeList[node2], 0, 0, 0, 0);
      AddAdjacency(nodeList[node2], nodeList[node1], 0, 0, 0, 0);
      paList[norm++]=node1;
      paList[norm++]=node2;
    }
  }

  /* Done */
  free(nodeList);
  free(paList);
  return root;
}

/*
  ---------------------------------------------------------------------
  Create an undirected block graph. Arguments:

  - ngroup: number of groups
  - gsize: a list with the group sizes
  - q: matrix of block-to-block connectivity probability
  - output_sw: 'v' for verbose
  - gen: random number generator

  ---------------------------------------------------------------------
*/
struct node_gra *
UndirectedBlockGraph(int ngroup, int *gsize, double **q,
		     char output_sw, gsl_rng *gen)
{
  int node1, node2;
  struct node_gra **nodeList = NULL;
  struct node_gra *root = NULL, *last = NULL;
  char label[MAX_LABEL_LENGTH];
  int i, S=0;

  /* Create header */
  last = root = CreateHeaderGraph();
  
  /* Create nodes */
  for (i=0; i<ngroup; i++)
    S += gsize[i];
  nodeList = (struct node_gra **)calloc(S, sizeof(struct node_gra *));
  S = 0;
  for (i=0; i<ngroup; i++) {
    for (node1=S; node1<S+gsize[i]; node1++) {
      sprintf(&label[0], "%d", node1+1);
      last = nodeList[node1] = CreateNodeGraph(last, &label[0]);
      last->inGroup = i;
    }
    S += gsize[i];
  }

  /* Create the links */
  for (node1=0; node1<S; node1++) {
    for (node2=node1+1; node2<S; node2++) {
      switch (output_sw) {
      case 'v':
	fprintf(stderr, "%s %s %g\n",
		nodeList[node1]->label, nodeList[node2]->label,
		q[(nodeList[node1])->inGroup][(nodeList[node2])->inGroup]);
	break;
      default:
	break;
      }
      if (gsl_rng_uniform(gen) <
	  q[(nodeList[node1])->inGroup][(nodeList[node2])->inGroup]) {
	AddAdjacency(nodeList[node1], nodeList[node2], 0, 0, 0, 0);
	AddAdjacency(nodeList[node2], nodeList[node1], 0, 0, 0, 0);
      }
    }
  }

  /* Done */
  free(nodeList);
  return root;
}

/*
  ---------------------------------------------------------------------
  Create a Girvan-Newman graph (PNAS, 2002) with community structure
  ---------------------------------------------------------------------
*/
struct node_gra *
GirvanNewmanGraph(int ngroup, int gsize, double kin, double kout,
		  char output_sw, gsl_rng *gen)
{
  double **q;
  int *gsizes;
  int g1, g2;
  struct node_gra *net=NULL;

  /* Prepare arrays necessary for UndirectedBlockGraph */
  q = allocate_d_mat(ngroup, ngroup);
  gsizes = allocate_i_vec(ngroup);
  for (g1=0; g1<ngroup; g1++) {
    gsizes[g1] = gsize;
    q[g1][g1] = (double)(kin) / (double)(gsize - 1);
    for (g2=g1+1; g2<ngroup; g2++) {
      q[g1][g2] = (double)(kout) / (double)(gsize * (ngroup - 1));
    }
  }

  /* Do it */
  net = UndirectedBlockGraph(ngroup, gsizes, q, output_sw, gen);

  /* Done */
  free_d_mat(q, ngroup);
  free_i_vec(gsizes);
  return net;
}
