/*
  main_betweenness.c
  $LastChangedDate: 2008-10-07 16:15:28 -0500 (Tue, 07 Oct 2008) $
  $Revision: 129 $
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include "graph.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL;
  gsl_rng *rand_gen;
  struct node_gra *net=NULL;
  int seed;

  /*
    ---------------------------------------------------------------------------
    Command line parameters
    ---------------------------------------------------------------------------
  */
  if (argc < 2) {
    printf("\nUse: netrandomize.out net_file seed\n\n");
    return -1;
  }
  netF = argv[1];
  seed = atoi(argv[2]);
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, seed);

  /*
    ---------------------------------------------------------------------------
    Build the network
    ---------------------------------------------------------------------------
  */
  infile = fopen(netF, "r");
  net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);

  /*
    ---------------------------------------------------------------------------
    Randomize the network
    ---------------------------------------------------------------------------
  */
  net = RandomizeSymmetricNetwork(net, 100, rand_gen);

  /*
    ---------------------------------------------------------------------------
    Output results and finish
    ---------------------------------------------------------------------------
  */
  FPrintNetAdjacencyList(stdout, net, 0, 1);
  RemoveGraph(net);
  gsl_rng_free(rand_gen);
  return 0;
}
