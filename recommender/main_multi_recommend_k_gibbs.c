/*
  -----------------------------------------------------------------------------
  Use:

  multi_recommend_k.out ratings_file_name query_file_name seed

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
#include "recommend.h"

int
main(int argc, char **argv)
{
  char *obsFile, *queFile;
  FILE *infile=NULL, *outfile=NULL;
  struct binet *ratings=NULL;
  gsl_rng *rand_gen;
  struct query **queries=NULL;
  int nquery, q;
  int seed;
  double **scores;
  int K, k;

  /* Command line parameters */
  if (argc < 4) {
    printf("\nUse: multi_recommend_k_gibbs K observation_file query_file seed\n\n");
    return -1;
  }
  K = atoi(argv[1]);
  obsFile = argv[2];
  queFile = argv[3];
  seed = atoi(argv[4]);
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, seed);

  /* Build the observations */
  infile = fopen(obsFile, "r");
  ratings = ReadRecommenderObservations(infile);
  fclose(infile);

  /* Get the number of queries end build the query set */
  nquery = CountLinesInFile(queFile);
  fprintf(stderr, ">> %d queries\n", nquery);
  infile = fopen(queFile, "r");
  queries = ReadQueries(infile, nquery, ratings);
  fclose(infile);

  /* Get the scores of the queries and print them */
  scores = GibbsMultiLinkScoreKState(K, ratings,
				     queries, nquery,
				     10000, rand_gen, 'v');
  fprintf(stdout, "\n\n>>> RESULTS\n\n");
  for (q=0; q<nquery; q++) {
    fprintf(stdout, "%s %s",
	    ((queries[q])->n1)->label,
	    ((queries[q])->n2)->label);
    for (k=0; k<K; k++)
      fprintf(stdout, " %lf", scores[k][q]);
    fprintf(stdout, "\n");
  }

  /* Finish */
  RemoveBipart(ratings);
  /* NEED TO FREE THE WHOLE QUERYSET!! */
/*   FreeQuery(the_query); */
  gsl_rng_free(rand_gen);
  return 0;
}
