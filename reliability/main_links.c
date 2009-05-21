/*
  main_links.c
  $LastChangedDate: 2008-03-25 18:02:36 -0500 (Tue, 25 Mar 2008) $
  $Revision: 124 $
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"
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
  struct prng *rand_gen;
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
    Get link reliabilities
    ---------------------------------------------------------------------------
  */
  newA = LinkScore(net, 0.0, 10000, rand_gen, 'v');
  printf("Done!\n");

  /*
    ---------------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------------
  */
  outfile1 = fopen("missing.dat", "w");
  outfile2 = fopen("bogus.dat", "w");
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
  printf("Done!\n");

  /*
    ---------------------------------------------------------------------------
    Finish
    ---------------------------------------------------------------------------
  */
  RemoveGraph(net);
  return 0;
}
