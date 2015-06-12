/*
  -----------------------------------------------------------------------------
  Use:

  multi_recommend_2.out ratings_file_name query_file_name seed

  seed is an integer, used to get the random number generator started
  -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "bipartite.h"
#include "conflict.h"

int
main(int argc, char **argv)
{
  char *obsFile, *queFile;
  FILE *infile=NULL, *outfile=NULL;
  struct binet *ratings=NULL;
  gsl_rng *rand_gen;
  int seed;
  int q, nquery;
  struct query **querySet;

  /* Command line parameters */
  if (argc < 3) {
    printf("\nUse: multi_recommend_2.out observation_file seed\n\n");
    return -1;
  }
  obsFile = argv[1];
  seed = atoi(argv[2]);
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, seed);

  /* Build the observation network */
  infile = fopen(obsFile, "r");
  ratings = ReadRecommenderObservations(infile);
  fclose(infile);
  nquery = CountNodes(ratings->net1) * CountNodes(ratings->net2);

  /* Get the scores and print them */
  querySet = AllLinkScore2State(ratings,
				10000, rand_gen, 'v', -1);
  outfile = fopen("link_scores.dat", "w");
  for (q=0; q<nquery; q++) {
    if (strcmp(querySet[q]->n1->label, querySet[q]->n2->label) != 0) {
      fprintf(outfile, "%s %s %lf\n",
	      (querySet[q]->n1)->label,
	      (querySet[q]->n2)->label,
	      querySet[q]->score);
    }
  }
  fclose(outfile);

  /* Finish */
  for (q=0; q<nquery; q++)
    FreeQuery(querySet[q]);
  free(querySet);
  RemoveBipart(ratings);
  gsl_rng_free(rand_gen);
  return 0;
}
