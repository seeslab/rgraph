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
  struct group *part2 = NULL;
  int S;

  while (1) {
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
    // Partition
    // ------------------------------------------------------------
    inFile = fopen("testpart.dat", "r");
    part = FCreatePartition(inFile);
    fclose(inFile);
    part = ClustersPartition(net);
    MapPartToNet(part, net);
/*     part = CreateEquiNPartitionSoft(4, 5); */
    FPrintPartition(stdout, part, 0);
    part2 = CreatePartitionFromInGroup(net);
    MapPartToNet(part2, net);
    FPrintPartition(stdout, part2, 0);

    // ------------------------------------------------------------
    // Free memory
    // ------------------------------------------------------------
    RemovePartition(part);
    RemoveGraph(net);
  }

  return 0;
}
