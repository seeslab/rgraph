#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <mpi.h>

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
    printf("\nUse: bipartmod_w_cl net_file_name seed iteration_factor cooling_factor modules_column\n\n");
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
  
  /* Initialize MPI */
  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

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

  if (world_rank == 0) {
    outF = fopen("modules_bipart_w_cl.dat", "w");
    FPrintPartition(outF, part, 0);
    fclose(outF);
  }

  /* Finalize MPI */
  MPI_Finalize();

  // Free memory
  // ------------------------------------------------------------
  RemovePartition(part);
  RemoveBipart(binet);
  gsl_rng_free(randGen);
}
