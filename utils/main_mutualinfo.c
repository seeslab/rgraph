#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include "graph.h"
#include "modules.h"

int main(int argc,char* argv[])
{
  FILE *inf = NULL;
  struct group *part1 = NULL;
  struct group *part2 = NULL;

  /*
    ------------------------------------------------------------
    Build the network and the partition
    ------------------------------------------------------------
  */
  //fprintf(stderr, "Creating the partition...\n");
  inf = fopen(argv[1], "r");
  part1 = FCreatePartition(inf);
  fclose(inf);

  //fprintf(stderr, "Creating the partition...\n");
  inf = fopen(argv[2], "r");
  part2 = FCreatePartition(inf);
  fclose(inf);

  /*
    ------------------------------------------------------------
    Calculate and print out the modularity
    ------------------------------------------------------------
  */

  printf("MI = %f\n", MutualInformation(part1, part2, 1));
      
  /*
    ------------------------------------------------------------
    Free memory
    ------------------------------------------------------------
  */
  RemovePartition(part1);
  RemovePartition(part2);

  return 0;
}
