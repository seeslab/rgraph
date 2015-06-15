#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include "tools.h"
#include "graph.h"
#include "layout.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *inFile;
  int i, nnet = argc - 1;
  struct node_gra *net[nnet];
  struct node_gra *net_sum = NULL, *new_net_sum = NULL;
  gsl_rng *rand_gen;

  /*
    ---------------------------------------------------------------------------
    Command line parameters
    ---------------------------------------------------------------------------
  */
  if (argc < 3) {
    printf("\nUse: multinetlayout.out net_file_1 net_file_2 [net_file_3...]\n\n");
    return -1;
  }

  /*
    ---------------------------------------------------------------------------
    Initialize the random number generator
    ---------------------------------------------------------------------------
  */
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, 1111);

  /*
    ---------------------------------------------------------------------------
    Build the networks
    ---------------------------------------------------------------------------
  */
  for (i=0; i<nnet; i++) {
    netF = argv[i+1];
    inFile = fopen(netF, "r");
    net[i] = FBuildNetwork(inFile, 0, 0, 0, 1);
    fclose(inFile);
  }

  /*
    ---------------------------------------------------------------------------
    Build and layout the sum network, and print the node coordinates
    ---------------------------------------------------------------------------
  */
  net_sum = AddTwoNetworks(net[0], net[1]);
  for (i=2; i<nnet; i++) {
    new_net_sum = AddTwoNetworks(net_sum, net[i]);
    RemoveGraph(net_sum);
    net_sum = new_net_sum;
  }
  MDGraphLayout(net_sum, 0.05, 0.05, 100000, rand_gen, 0);
  PrintNodeCoordinates(stdout, net_sum);

  /*
    ---------------------------------------------------------------------------
    Free memory
    ---------------------------------------------------------------------------
  */
  for (i=0; i<nnet; i++)
    RemoveGraph(net[i]);
  RemoveGraph(net_sum);
  gsl_rng_free(rand_gen);
  return 0;
}
