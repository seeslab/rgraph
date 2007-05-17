#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "graph.h"

main()
{
  FILE *inFile;
  struct node_gra *net = NULL;
  int S;

  // ------------------------------------------------------------
  // Build the network
  // ------------------------------------------------------------
  inFile = fopen("test.dat", "r");
  net = FBuildNetwork(inFile, 0, 1, 0, 1);
  fclose(inFile);

  S = CountNodes(net);
  printf("The net has %d nodes\n", S);

  RemoveGraph(net);
}
