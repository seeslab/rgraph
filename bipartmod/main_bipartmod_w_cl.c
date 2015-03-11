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
  int to_file = 0;
  int from_file = 1;
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


  file_name = argv[1];
  seed = atoi(argv[2]);
  fac = atof(argv[3]);
  Ts = atof(argv[4]);
  invert = atoi(argv[5]);
  to_file = atoi(argv[6]);
  from_file = atoi(argv[7]);
  
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
  if (from_file == 1) {
	inF = fopen(file_name, "r");
	binet = FBuildNetworkBipart(inF, 1, 0);
	fclose(inF);
  }
  else
	binet = FBuildNetworkBipart(stdin, 1, 0);

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
  if (to_file == 1)	{
	outF = fopen("modules_bipart_w_cl.dat", "w");
    FPrintPartition(outF, part, 0);
	fclose(outF);
  }
  else{
	printf("### Results:\n");
	FPrintPartition(stdout, part, 0);
  }
  // Free memory
  // ------------------------------------------------------------
  RemovePartition(part);
  RemoveBipart(binet);
  gsl_rng_free(randGen);
}
