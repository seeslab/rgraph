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
/*     int S; */
/*     S = CountNodes(net); */
/*     printf("The net has %d nodes\n", S); */
/*     printf("C = %g\n", ClusteringCoefficient(net)); */
/*     FPrintPajekFile("test.net", net, 0, 0, 1); */

    // ------------------------------------------------------------
    // Partition
    // ------------------------------------------------------------
/*     struct group *part = NULL; */
/*     inFile = fopen("testpart.dat", "r"); */
/*     part = FCreatePartition(inFile); */
/*     fclose(inFile); */
/*     MapPartToNet(part, net); */
/*     fprintf(stderr, "\n"); */
/*     FPrintPartition(stdout, part, 0); */
/*     RemovePartition(part); */

/*     struct node_gra *net2; */
/*     net2 = BuildNetFromPart(part->next); */
/*     RemoveGraph(net2); */

/*     struct group *part2 = NULL; */
/*     part2 = ClustersPartition(net); */
/*     MapPartToNet(part2, net); */
/*     FPrintPartition(stdout, part2, 0); */
/*     RemovePartition(part2); */

    struct group *part2 = NULL;
    part2 = SACommunityIdent(net, 2. / CountNodes(net), 0.0, 0.995,
			     2.0, 2, 'o', 1, 'v', randGen);
    MapPartToNet(part2, net);
    FPrintPartition(stdout, part2, 0);
    RemovePartition(part2);

    // ------------------------------------------------------------
    // Free memory
    // ------------------------------------------------------------
    RemoveGraph(net);
    prng_free(randGen);
  }

  return 0;
}
