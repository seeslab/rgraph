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
  fprintf(stderr, "Creating the network...\n");
  binet = BuildModularBipartiteNetwork(NULL, 32, 4, NULL, 64,
				       7, 7, -1.,
				       0.85, randGen);

  fprintf(stderr, "Inverting the network twice...\n");
  InvertBinet(binet);
  InvertBinet(binet);

  fprintf(stderr, "Randomizing the network...\n");
  RandomizeBinet(binet, 100, randGen);

  fprintf(stderr, "Projecting the network...\n");
  struct node_gra *projection;
  projection = ProjectBinet(binet);
  RemoveGraph(projection);

  FPrintPajekFileBipart ("binet.net", binet, 0, 0);

  /*
    ------------------------------------------------------------
    Network properties
    ------------------------------------------------------------
  */
  printf("L = %d\n", NLinksBinet(binet));
      
  /*
    ------------------------------------------------------------
    Free memory
    ------------------------------------------------------------
  */
  RemoveBinet(binet);
  prng_free(randGen);

  return 0;
}
