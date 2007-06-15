#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"

#include "graph.h"
#include "models.h"
#include "modules.h"

int main()
{
  struct node_gra *net = NULL;
  struct prng *randGen;
  int repeat = 1;

  // --------------------------------------------------------------
  // Initialize the random number generator
  // --------------------------------------------------------------
    randGen = prng_new("mt19937(1111)");
    prng_seed(randGen, 3333);

  while (repeat) {

    // ------------------------------------------------------------
    // Build the network
    // ------------------------------------------------------------
/*     FILE *inf = fopen("karate.dat", "r"); */
/*     net = FBuildNetwork(inf, 0, 1, 0, 1); */
/*     fclose(inf); */

/*     net = FBuildNetwork(stdin, 0, 1, 0, 1); */

    net = ERGraph(50, 0.06, randGen);

    // ------------------------------------------------------------
    // Network properties
    // ------------------------------------------------------------
    printf("S = %d\n", CountNodes(net));
    printf("C = %g\n", ClusteringCoefficient(net));
    printf("A = %g\n", Assortativity(net));
    printf("P = %g\n", AverageInverseDistance(net));

/*     FPrintPajekFile("test.net", net, 0, 0, 1); */

    // ------------------------------------------------------------
    // Partition
    // ------------------------------------------------------------
/*     struct group *part = NULL; */
/*     part = SACommunityIdent(net, 2. / CountNodes(net), 0.0, 0.995, */
/* 			     2.0, 2, 'o', 1, 'n', randGen); */
/*     MapPartToNet(part, net); */
/* /\*     FPrintPartition(stdout, part, 0); *\/ */
/*     printf("M = %g\n", Modularity(part)); */
/*     RemovePartition(part); */

    // ------------------------------------------------------------
    // Free memory
    // ------------------------------------------------------------
    RemoveGraph(net);

/*     repeat = 0; */
  }

  prng_free(randGen);

  return 0;
}
