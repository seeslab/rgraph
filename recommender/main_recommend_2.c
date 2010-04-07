/*
  -----------------------------------------------------------------------------
  Use:

  state2.out network_file_name node1 node2 seed

  seed is an integer, used to get the random number generator started
  -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"
#include "graph.h"
#include "recommend.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL, *outfile=NULL;
  struct binet *binet = NULL;
  struct prng *rand_gen;
  int seed;

  /* Command line parameters */
  if (argc < 3) {
    printf("\nUse: recommend_2.out net_file seed\n\n");
    return;
  }
  netF = argv[1];
  seed = atoi(argv[2]);

  /* Build the network */
  infile = fopen(netF, "r");
  binet = FBuildNetworkBipart(infile, 0, 0);
  fclose(infile);

  /* Get the score of the query */
  rand_gen = prng_new("mt19937(1111)");
  prng_seed(rand_gen, seed);
  LinkScore2State(binet, the_query, 10000, rand_gen, 'd');

  /* Finish */
  RemoveBipart(net);

  return 0;
}
