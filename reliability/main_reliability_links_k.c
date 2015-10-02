/*
  main_links.c
  $LastChangedDate: 2008-03-25 18:02:36 -0500 (Tue, 25 Mar 2008) $
  $Revision: 124 $
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "missing.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL;
  FILE *outfile1=NULL, *outfile2=NULL;
  struct node_gra *net=NULL;
  gsl_rng *rand_gen;
  double **scores;
  struct node_gra *p1, *p2;
  long int seed;
  int K, k;
  struct query **queries=NULL;
  int nquery=0, q;

  /*
    ---------------------------------------------------------------------------
    Command line parameters
    ---------------------------------------------------------------------------
  */
  if (argc < 3) {
    printf("\nUse: reliability_links_k K net_file seed\n\n");
    return -1;
  }
  K = atoi(argv[1]);
  netF = argv[2];
  seed = atoi(argv[3]);
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, seed);

  /*
    ---------------------------------------------------------------------------
    Build the network
    ---------------------------------------------------------------------------
  */
  infile = fopen(netF, "r");
  net = FBuildNetwork(infile, 1, 0, 0, 1);
  fclose(infile);

  /*
    ---------------------------------------------------------------------------
    Get link reliabilities
    ---------------------------------------------------------------------------
  */
  scores = LSMultiLinkScoreKState(K, net, 10000, rand_gen, 'v', -1,
				  &queries, &nquery);

  /*
    ---------------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------------
  */
  outfile1 = fopen("unobserved_k.dat", "w");
  outfile2 = fopen("observed_k.dat", "w");
  for (q=0; q<nquery; q++) {
    if (IsThereLink(queries[q]->n1, queries[q]->n2) == 0) {
      fprintf(outfile1, "%s %s",
	      ((queries[q])->n1)->label,
	      ((queries[q])->n2)->label);
      for (k=0; k<K; k++)
	fprintf(outfile1, " %lf", scores[k][q]);
      fprintf(outfile1, "\n");
    }
    else {
      fprintf(outfile2, "%s %s",
	      ((queries[q])->n1)->label,
	      ((queries[q])->n2)->label);
      for (k=0; k<K; k++)
	fprintf(outfile2, " %lf", scores[k][q]);
      fprintf(outfile2, "\n");
    }
  }
  fclose(outfile1);
  fclose(outfile2);

  /*
    ---------------------------------------------------------------------------
    Finish
    ---------------------------------------------------------------------------
  */
  RemoveGraph(net);
  gsl_rng_free(rand_gen);
  return 0;
}
