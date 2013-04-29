#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include "bipartite.h"
#include "graph.h"
#include "modules.h"

int main(int argc,char* argv[])
{
  FILE *inf = NULL;
  struct binet *binet = NULL;
  struct group *part = NULL;

  /*
    ------------------------------------------------------------
    Build the network and the partition
    ------------------------------------------------------------
  */
  //fprintf(stderr, "Creating the network...");
  inf = fopen(argv[1], "r");
  binet = FBuildNetworkBipart(inf, 0, 0);
  fclose(inf);
  //fprintf(stderr, "done.\n");

  //fprintf(stderr, "Creating the partition...");
  inf = fopen(argv[2], "r");
  part = FCreatePartition(inf);
  fclose(inf);
  //fprintf(stderr, "done.\n");

  //fprintf(stderr, "Inverting the network...");
  if(atoi(argv[3])==1){
    InvertBipart(binet);
  }
  //fprintf(stderr, "done.\n");

  //fprintf(stderr, "Mapping the partition...");
  MapPartToNet(part, binet->net1);
  //fprintf(stderr, "done.\n");

  /*
    ------------------------------------------------------------
    Calculate and print out the modularity
    ------------------------------------------------------------
  */

  //fprintf(stderr, "Calculating the weighted bipartite modularity...");
  double M = ModularityBipart(binet, part);
  //fprintf(stderr, "done.\n");

  printf("M = %f\n", M);
      
  /*
    ------------------------------------------------------------
    Free memory
    ------------------------------------------------------------
  */
  RemovePartition(part);
  RemoveBipart(binet);

  return 0;
}
