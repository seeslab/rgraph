#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"

#include "graph.h"
#include "modules.h"

int main()
{
  FILE *inFile;
  struct node_gra *net = NULL;
  struct group *part = NULL;
  struct group *part2 = NULL;
  int S;
  struct prng *randGen;

  while (1) {

    // ------------------------------------------------------------
    // Initialize the random number generator
    // ------------------------------------------------------------
    randGen = prng_new("mt19937(1111)");
    prng_seed(randGen, 1111);

    // ------------------------------------------------------------
    // Build the network
    // ------------------------------------------------------------
    inFile = fopen("test.dat", "r");
    net = FBuildNetwork(inFile, 0, 1, 0, 1);
    fclose(inFile);

    // ------------------------------------------------------------
    // Network
    // ------------------------------------------------------------
/*     S = CountNodes(net); */
/*     printf("The net has %d nodes\n", S); */
/*     printf("C = %g\n", ClusteringCoefficient(net)); */
/*     FPrintPajekFile("test.net", net, 0, 0, 1); */

    // ------------------------------------------------------------
    // Partition
    // ------------------------------------------------------------
/*     inFile = fopen("testpart.dat", "r"); */
/*     part = FCreatePartition(inFile); */
/*     fclose(inFile); */
/*     part = ClustersPartition(net); */
/*     MapPartToNet(part, net); */
/*     part = CreateEquiNPartitionSoft(4, 5); */
/*     FPrintPartition(stdout, part, 0); */

/*     part2 = CreatePartitionFromInGroup(net); */

    part2 = SACommunityIdent(net, 2. / CountNodes(net), 0.0, 0.995,
			     2.0, 2, 'o', 1, 'v', randGen);
    MapPartToNet(part2, net);
    FPrintPartition(stdout, part2, 0);

    // ------------------------------------------------------------
    // Free memory
    // ------------------------------------------------------------
/*     RemovePartition(part); */
    RemovePartition(part2);
    RemoveGraph(net);
    prng_free(randGen);
  }

  return 0;
}
