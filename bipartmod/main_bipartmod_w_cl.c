#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>
#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "bipartite.h"

main(int argc, char **argv)
{
  FILE *outF, *inF;
  int rgm;
  int seed = 1111;
  struct binet *binet = NULL;
  struct group *part = NULL;
  gsl_rng *randGen;
  double Ti, Tf, Ts, fac;
  char *file_name;
  int invert;

  /*
    ------------------------------------------------------------
    Prompt for user-defined parameters
    ------------------------------------------------------------
  */
  if (argc < 6) {
    printf("\nUsage: bipartmod_cl net_file_name seed iteration_factor cooling_factor modules_column\n\n"
		   "\t net_file_name: Name of the network file\n"
		   "\t seed: Random number seed (POSITIVE Integer) \n"
		   "\t iteraction_factor: Iteration factor (recommended 1.0)\n"
		   "\t cooling_factor: Cooling factor (recommended 0.950-0.995)\n "
		   "\t module_column: Find modules for the first column (0) or second columnd (1)\n\n");
    return -1;
  }

  //printf("\n# Enter the name of the network file: ");
  file_name = argv[1];
  
  //printf("\n# Enter random number seed (POSITIVE integer): ");
  seed = atoi(argv[2]);
  
  //printf("\n# Enter iteration factor (recommended 1.0): ");
  fac = atof(argv[3]);
  
  //printf("\n# Enter the cooling factor (recommended 0.950-0.995): ");
  Ts = atof(argv[4]);
  
  //printf("\n# Find modules from first column (0) or second columnd (1): ");
  invert = atoi(argv[5]);
  
  /*
    ------------------------------------------------------------------
    Initialize the random number generator
    ------------------------------------------------------------------
  */
  randGen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(randGen, seed);

  /*
    ------------------------------------------------------------------
    Build the network
    ------------------------------------------------------------------
  */
  inF = fopen(file_name, "r");
  binet = FBuildNetworkBipart(inF, 1, 0);
  fclose(inF);
  if (invert == 1)
    InvertBipart(binet);

  /*
    ------------------------------------------------------------------
    Find the modules using the bipartite network
    ------------------------------------------------------------------
  */
  Ti = 1. / (double)CountNodes(binet->net1);
  Tf = 0.;

  part = SACommunityIdentBipartWeighted(binet,
				Ti, Tf, Ts, fac,
				0, 'o', 1, 'm',
				randGen);
  
  outF = fopen("modules_bipart_w_cl.dat", "w");
  FPrintPartition(outF, part, 0);
  fclose(outF);

  // Free memory
  // ------------------------------------------------------------
  RemovePartition(part);
  RemoveBipart(binet);
  gsl_rng_free(randGen);
}
