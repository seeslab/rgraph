/*
  models.c
  $LastChangedDate$
  $Revision$
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "prng.h"

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
ERGraph(int S, double p, struct prng *gen)
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
      if (prng_get_next(gen) < p) {
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
PAGraph(int S, int m, struct prng *gen)
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
	node2 = paList[(int)(prng_get_next(gen) * norm)];
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
