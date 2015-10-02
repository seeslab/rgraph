/*
  main_reliability_links_mb.c
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "multiblock.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL;
  FILE *outfileAND=NULL, *outfileOR=NULL;
  struct node_gra *net=NULL;
  gsl_rng *rand_gen;
  double **newA_AND;
  struct node_gra *p1, *p2;
  long int seed;
  char outFileNameOR[200], outFileNameAND[200];

  /*
    ---------------------------------------------------------------------------
    Command line parameters
    ---------------------------------------------------------------------------
  */
  if (argc < 2) {
    printf("\nUse: links net_file seed\n\n");
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
    Get link reliabilities
    ---------------------------------------------------------------------------
  */
  newA_AND = LinkScoreMB(net, 0.0, 10000, rand_gen, 'q');

  /*
    ---------------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------------
  */
  strcpy(outFileNameAND, netF);
  strcat(outFileNameAND, ".AND_scores");
  outfileAND = fopen(outFileNameAND, "w");
  p1 = net;
  while ((p1 = p1->next) != NULL) {
    p2 = p1;
    while ((p2 = p2->next) != NULL) {
      fprintf(outfileAND,
	      "%g %s %s\n", newA_AND[p1->num][p2->num], p1->label, p2->label);
    }
  }
  fclose(outfileAND);

  /*
    ---------------------------------------------------------------------------
    Finish
    ---------------------------------------------------------------------------
  */
  RemoveGraph(net);
  gsl_rng_free(rand_gen);
  return 0;
}
