/*
  -----------------------------------------------------------------------------
  Network reconstruction based on the stochastic blockmodel sampling approach

  Use:

  reconstruct.out network_file_name seed

  seed is an integer, used to get the random number generator started
  -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"
#include "graph.h"
#include "missing.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL, *outfile=NULL;
  struct node_gra *net=NULL, *netRec=NULL;
  struct prng *rand_gen;
  int seed;

  /* Command line parameters */
  if (argc < 3) {
    printf("\nUse: reconstruct.out net_file seed\n\n");
    return;
  }
  netF = argv[1];
  seed = atoi(argv[2]);

  /* Build the network */
  infile = fopen(netF, "r");
  net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);

  /* Reconstruct the network */
  rand_gen = prng_new("mt19937(1111)");
  prng_seed(rand_gen, seed);
  netRec = NetReconstruct(net, rand_gen);

  /* Print the reconstruction */
  outfile = fopen("net_reconstructed.dat", "w");
  FPrintNetAdjacencyList(outfile, netRec, 0, 1);
  fclose(outfile);

  /* Finish */
  RemoveGraph(net);
  RemoveGraph(netRec);

  return 0;
}
