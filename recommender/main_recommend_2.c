/*
  -----------------------------------------------------------------------------
  Use:

  recommend_2.out network_file_name node1 node2 seed

  seed is an integer, used to get the random number generator started
  -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"
#include "graph.h"
#include "bipartite.h"
#include "recommend.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL, *outfile=NULL;
  struct binet *binet=NULL;
  struct prng *rand_gen;
  struct query *the_query=NULL;
  void *dict1, *dict2;
  char *node1, *node2;
  int seed;
  double score;

  /* Command line parameters */
  if (argc < 3) {
    printf("\nUse: recommend_2.out net_file node1 node2 seed\n\n");
    return;
  }
  netF = argv[1];
  node1 = argv[2];
  node2 = argv[3];
  seed = atoi(argv[4]);
  rand_gen = prng_new("mt19937(1111)");
  prng_seed(rand_gen, seed);

  /* Build the network */
  infile = fopen(netF, "r");
  binet = FBuildNetworkBipart(infile, 0, 0);
  fclose(infile);

  /* Get the query */
  dict1 = MakeLabelDict(binet->net1);
  dict2 = MakeLabelDict(binet->net2);
  the_query = CreateQuery(GetNodeDict(node1, dict1),
			  GetNodeDict(node2, dict2));

  /* Get the score of the query */
  score = LinkScore2State(binet, the_query, 10000, rand_gen, 'd');
  fprintf(stdout, ">>> SCORE: %g\n", score);

  /* Finish */
  RemoveBipart(binet);
  FreeQuery(the_query);
  return 0;
}
