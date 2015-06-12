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
  double **newA;
  struct node_gra *p1, *p2;
  long int seed;

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
  newA = LinkScore(net, 0.0, 10000, rand_gen, 'v');

  /*
    ---------------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------------
  */
  outfile1 = fopen("missing.dat", "w");
  outfile2 = fopen("spurious.dat", "w");
  p1 = net;
  while ((p1 = p1->next) != NULL) {
    p2 = p1;
    while ((p2 = p2->next) != NULL) {
      if (IsThereLink(p1, p2) == 0) {
	fprintf(outfile1,
		"%g %s %s\n", newA[p1->num][p2->num], p1->label, p2->label);
      }
      else {
	fprintf(outfile2,
		"%g %s %s\n", newA[p1->num][p2->num], p1->label, p2->label);
      }
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
