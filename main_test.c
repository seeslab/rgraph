#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"

#include "graph.h"
#include "modules.h"

int main()
{
  struct node_gra *net = NULL;
  struct prng *randGen;
  int repeat = 1;

  while (repeat) {

    // ------------------------------------------------------------
    // Initialize the random number generator
    // ------------------------------------------------------------
    randGen = prng_new("mt19937(1111)");
    prng_seed(randGen, 1111);

    // ------------------------------------------------------------
    // Build the network
    // ------------------------------------------------------------
/*     FILE *inFile; */
/*     inFile = fopen("test.dat", "r"); */
/*     net = FBuildNetwork(inFile, 0, 1, 0, 1); */
/*     fclose(inFile); */

    net = FBuildNetwork(stdin, 0, 1, 0, 1);

    // ------------------------------------------------------------
    // Network
    // ------------------------------------------------------------
/*     int S; */
/*     S = CountNodes(net); */
/*     printf("The net has %d nodes\n", S); */
/*     printf("C = %g\n", ClusteringCoefficient(net)); */
/*     FPrintPajekFile("test.net", net, 0, 0, 1); */

    // ------------------------------------------------------------
    // Partition
    // ------------------------------------------------------------
    struct group *part = NULL;
    part = SACommunityIdent(net, 2. / CountNodes(net), 0.0, 0.995,
			     2.0, 2, 'o', 1, 'n', randGen);
    MapPartToNet(part, net);
    FPrintPartition(stdout, part, 0);
    RemovePartition(part);

    // ------------------------------------------------------------
    // Free memory
    // ------------------------------------------------------------
    RemoveGraph(net);
    prng_free(randGen);

    repeat = 0;
  }

  return 0;
}
