/*
  main_missing.c
  $LastChangedDate: 2008-03-25 18:02:36 -0500 (Tue, 25 Mar 2008) $
  $Revision: 124 $
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"
#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "missing.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL;
  struct node_gra *net=NULL;
  struct prng *rand_gen;
  double **newA;
  struct node_gra *p1, *p2;

  /*
    ---------------------------------------------------------------------------
    Command line parameters
    ---------------------------------------------------------------------------
  */
  if (argc < 2) {
    printf("\nUse: missing.out net_file\n\n");
    return;
  }
  netF = argv[1];
  rand_gen = prng_new("mt19937(1111)");
  prng_seed(rand_gen, 3333);

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
    Missing
    ---------------------------------------------------------------------------
  */
  newA = LinkScore(net, 0.0, 10000, rand_gen);

  p1 = p2 = net;
  while ((p1 = p1->next) != NULL) {
    p2 = p1;
    while ((p2 = p2->next) != NULL) {
      printf("\"%s-%s\" %g\n", p1->label, p2->label, newA[p1->num][p2->num]);
    }
  }

  /*
    ---------------------------------------------------------------------------
    Finish
    ---------------------------------------------------------------------------
  */
  RemoveGraph(net);
  return 0;
}
