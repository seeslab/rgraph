#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "prng.h"

#include "graph.h"
#include "models.h"

// ---------------------------------------------------------------------
// Create an Erdos-Renyi random graph
// ---------------------------------------------------------------------
struct node_gra *
ERGraph(int S, double p, struct prng *gen)
{
  int node1, node2;
  struct node_gra **nodeList = NULL;
  struct node_gra *root = NULL, *last = NULL;
  char label[MAX_LABEL_LENGTH];

  // Create header
  last = root = CreateHeaderGraph();
  
  // Create nodes
  nodeList = (struct node_gra **)calloc(S, sizeof(struct node_gra *));
  for (node1=0; node1<S; node1++) {
    sprintf(&label[0], "%d", node1+1);
    last = nodeList[node1] = CreateNodeGraph(last, &label[0]);
  }

  // Create the links
  for (node1=0; node1<S; node1++) {
    for (node2=node1+1; node2<S; node2++) {
      if (prng_get_next(gen) < p) {
	AddAdjacency(nodeList[node1], nodeList[node2], 0, 0, 0, 0);
	AddAdjacency(nodeList[node2], nodeList[node1], 0, 0, 0, 0);
      }
    }
  }

  // Done
  free(nodeList);
  return root;
}
