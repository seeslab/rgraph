#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "prng.h"

#include "graph.h"

// ---------------------------------------------------------------------
// Create an Erdos-Renyi random graph
// ---------------------------------------------------------------------
struct node_gra *ERGraph(int S, double p, struct prng *gen)
{
  int node1, node2;
  struct node_gra **list;
  struct node_gra *root = NULL, *last = NULL;
  char label[MAX_LABEL_LENGTH];

  // Create header
  last = root = CreateHeaderGraph();
  
  // Create nodes
  for (node1=0; node1<S; node1++)
    sprintf(&label[0], "%d", node1+1);
    last = list[node1] = CreateNodeGraph(last, &label[0]);

  // Create the links
  for (node1=0; node1<S; node1++) {
    for (node2=node1+1; node2<S; node2++) {
      if (prng_get_next(gen) < p) {
	AddAdjacency(list[node1], list[node2], 0, 0, 0, 0);
	AddAdjacency(list[node2], list[node1], 0, 0, 0, 0);
      }
    }
  }

  // Done
  return root;
}
