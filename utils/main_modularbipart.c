/*
  main_modularbipart.c
  $LastChangedDate: 2007-06-29 14:58:36 -0500 (Fri, 29 Jun 2007) $
  $Revision: 70 $
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include "graph.h"
#include "bipartite.h"

int
main(int argc, char **argv)
{
  struct binet *binet = NULL;
  gsl_rng *randGen;
  int nMod, modSize, nTeam, teamSize, seed;
  double teamHomog;

  /* Command line parameters */
  nMod = atoi(argv[1]);
  modSize = atoi(argv[2]);
  nTeam = atoi(argv[3]);
  teamSize = atoi(argv[4]);
  teamHomog = atof(argv[5]);
  seed = atoi(argv[6]);
  
  /* Initialize the random number generator */
  randGen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(randGen, seed);


  /* Build the network */
  binet = BuildModularBipartiteNetwork(NULL, modSize, nMod, NULL,
				       nTeam, modSize, modSize, -1.,
				       teamHomog, randGen);
  FPrintPajekFileBipart ("binet.net", binet, 0, 0);

  /* Free memory */
  RemoveBipart(binet);
  gsl_rng_free(randGen);

  return 0;
}
