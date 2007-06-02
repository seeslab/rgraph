#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "graph.h"
#include "modules.h"

int main()
{
  FILE *inFile;
  struct node_gra *net = NULL;
  struct group *part = NULL;

  int S;

  // ------------------------------------------------------------
  // Build the network
  // ------------------------------------------------------------
  inFile = fopen("test.dat", "r");
  net = FBuildNetwork(inFile, 0, 1, 0, 1);
  fclose(inFile);

  // ------------------------------------------------------------
  // Network
  // ------------------------------------------------------------
  S = CountNodes(net);
  printf("The net has %d nodes\n", S);
  printf("C = %g\n", ClusteringCoefficient(net));
  FPrintPajekFile("test.net", net, 0, 0, 1);

  // ------------------------------------------------------------
  // Partitions
  // ------------------------------------------------------------
  inFile = fopen("testpart.dat", "r");
  part = FBuildPartition(inFile);
  fclose(inFile);
  printf("First node: %s\n", part->next->nodeList->next->nodeLabel);
  MapPartToNet(part, net);
  FPrintPartition(stdout, part, 0);

  // ------------------------------------------------------------
  // Free memory
  // ------------------------------------------------------------
  RemoveGraph(net);
  RemovePartition(part);

  return 0;
}
