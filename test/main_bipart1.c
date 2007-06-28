/*
  main_bipart1.c
  $LastChangedDate$
  $Revision$
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"

#include "graph.h"
#include "bipartite.h"

int main()
{
  struct binet *binet = NULL;
  struct prng *randGen;

  /*
    --------------------------------------------------------------
    Initialize the random number generator
    --------------------------------------------------------------
  */
  randGen = prng_new("mt19937(1111)");
  prng_seed(randGen, 3333);

  /*
    ------------------------------------------------------------
    Build the network
    ------------------------------------------------------------
  */
  binet = BuildModularBipartiteNetwork(int *mod_size,
				       int nodpermod,
				       int nmod,
				       double *col_prob,
				       int S2,
				       int mmin, int mmax,
				       double geom_p,
				       double p,
				       struct prng *gen);

  /*
    ------------------------------------------------------------
    Network properties
    ------------------------------------------------------------
  */
  printf("S = %d\n", CountNodes(net));
  printf("C = %g\n", ClusteringCoefficient(net));
  printf("A = %g\n", Assortativity(net));
  printf("P = %g\n", AverageInverseDistance(net));
      
  /*
    ------------------------------------------------------------
    Partition
    ------------------------------------------------------------
  */
  struct group *part = NULL;
  part = SACommunityIdent(net, 2. / CountNodes(net), 0.0, 0.995,
			  2.0, 2, 'o', 1, 'n', randGen);
  MapPartToNet(part, net);
  printf("M = %g\n", Modularity(part));
  RemovePartition(part);

  /*
    ------------------------------------------------------------
    Free memory
    ------------------------------------------------------------
  */
  RemoveGraph(net);
  prng_free(randGen);

  return 0;
}
