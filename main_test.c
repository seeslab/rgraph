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
  while (1) {
    inFile = fopen("testpart.dat", "r");
    part = FBuildPartition(inFile);
    fclose(inFile);
    MapPartToNet(part, net);
    FPrintPartition(stdout, part, 0);
    RemovePartition(part);
  }

  // ------------------------------------------------------------
  // Free memory
  // ------------------------------------------------------------
  RemoveGraph(net);

  return 0;
}
