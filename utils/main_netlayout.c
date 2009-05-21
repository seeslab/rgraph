#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "prng.h"
#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "missing.h"
#include "matrix.h"

int
main(int argc, char **argv)
{
  long seed;
  char *netF;
  int method;
  FILE *inFile;
  struct node_gra *net = NULL;
  struct prng *rand_gen;
  int S, L;
  int nsteps;
  double step, damp;

  /*
    ---------------------------------------------------------------------------
    Command line parameters
    ---------------------------------------------------------------------------
  */
  if (argc < 7) {
    printf("\nUse: netlayout.out net_file_name seed step(-1) nsteps damping(-1) 1(2d)/2(2dp)/3(3d)\n\n");
    return;
  }

  netF = argv[1];
  seed = atoi(argv[2]);
  step = atof(argv[3]);
  if (step < 0.0 )
    step = 0.05; // Default step
  nsteps = atoi(argv[4]);
  damp = atof(argv[5]);
  if (damp < 0.0 )
    damp = 0.05; // Default damping
  method = atoi(argv[6]);

  /*
    ---------------------------------------------------------------------------
    Initialize the random number generator
    ---------------------------------------------------------------------------
  */
  rand_gen = prng_new("mt19937(1111)");
  prng_seed(rand_gen, seed);

  /*
    ---------------------------------------------------------------------------
    Build the network
    ---------------------------------------------------------------------------
  */
  inFile=fopen(netF, "r");
  net = FBuildNetwork(inFile, 0, 0, 0, 1);
  fclose(inFile);
  S = CountNodes(net);
  L = TotalNLinks(net, 1);

  /*
    ---------------------------------------------------------------------------
    Layout the graph
    ---------------------------------------------------------------------------
  */
  if (method == 3) {
    MDGraphLayout3D(net, damp, step, nsteps, rand_gen, 0);
  }
  else if (method == 2) {
    MDGraphLayout2Dp(net, damp, step, nsteps, rand_gen, 0);
  }
  else {
    MDGraphLayout(net, damp, step, nsteps, rand_gen, 0);
  }

  /*
    ---------------------------------------------------------------------------
    Output coordinates
    ---------------------------------------------------------------------------
  */
  PrintNodeCoordinates(stdout, net);

  /*
    ---------------------------------------------------------------------------
    Free memory
    ---------------------------------------------------------------------------
  */
  RemoveGraph(net);
}
