/*
  main_betweenness.c
  $LastChangedDate: 2008-10-07 16:15:28 -0500 (Tue, 07 Oct 2008) $
  $Revision: 129 $
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "graph.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL;
  struct prng *rand_gen;
  struct node_gra *net=NULL;
  int seed;

  /*
    ---------------------------------------------------------------------------
    Command line parameters
    ---------------------------------------------------------------------------
  */
  if (argc < 2) {
    printf("\nUse: netrandomize.out net_file seed\n\n");
    return;
  }
  netF = argv[1];
  seed = atoi(argv[2]);
  rand_gen = prng_new("mt19937(1111)");
  prng_seed(rand_gen, seed);

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
  return 0;
}
