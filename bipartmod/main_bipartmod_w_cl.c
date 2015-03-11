#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

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
  int from_file = 0;
  struct binet *binet = NULL;
  struct group *part = NULL;
  gsl_rng *randGen;
  double Ti, Tf;
  double Ts = 0.97;
  double fac = 1.0;
  char *file_name;
  char *file_name_out;
  int invert = 0;
  int weighted = 0;

  extern char *optarg;
  extern int optind;
  int c; 
  
  /*
    ------------------------------------------------------------
    Arguments parsing
    ------------------------------------------------------------
  */
  
  if (argc == 1) {
    printf("\nUsage: bipartmod_cl [-f file] [-o file] [-s seed] [-i iter] [-c cool] [-pw]\n\n"
		   "\t -f file: Name of the input network file (default: standard input)\n"
		   "\t -o file: Name of the output file (default: standard output)\n"
		   "\t -s seed: Random number seed (POSITIVE Integer, default 1111) \n"
		   "\t -i iter: Iteration factor (recommended 1.0, default 1.0)\n"
		   "\t -c cool: Cooling factor (recommended 0.950-0.995, default 0.97)\n "
		   "\t -p : Find modules for the second column (default: first) \n"
		   "\t -w : Use the third column of the input file and the weighted version of the algorithm\n");
	
    return -1;
  }
  else{	
  while ((c = getopt(argc, argv, "wpf:s:i:c:o:")) != -1)
	switch (c) {
	case 'f':
	  from_file = 1;
	  file_name = optarg;
	  //printf("input set to %s \n", file_name);
	  break;
	case 's':
	  seed = atoi(optarg);
	  //printf("seed set to %d \n", seed);
	  break;
	case 'i':
	  fac = atof(optarg);
	  //printf("iteration factor set to %f \n", fac);
	  break;
	case 'c':
	  Ts = atof(optarg);
	  //printf("cooling factor set to %f. \n", Ts);
	  break;
	case 'p':
	  invert = 1;
	  //printf("Looking at second column.\n");
	  break;
	case 'o':
	  to_file = 1;
	  file_name_out = optarg;
	  //printf("Output file set to %s \n",file_name_out);
	  break;
	case 'w':
	  weighted = 1;
	  break;
	}
  }


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
	binet = FBuildNetworkBipart(inF, weighted, 0);
	fclose(inF);
  }
  else
	binet = FBuildNetworkBipart(stdin, weighted, 0);

  if (invert == 1)
    InvertBipart(binet);

  /*
    ------------------------------------------------------------------
    Find the modules using the bipartite network
    ------------------------------------------------------------------
  */
  Ti = 1. / (double)CountNodes(binet->net1);
  Tf = 0.;

  if (weighted == 1){
	part = SACommunityIdentBipartWeighted(binet,
										  Ti, Tf, Ts, fac,
										  0, 'o', 1, 'm',
										  randGen);
  }
  else{
	part = SACommunityIdentBipart(binet,
								  Ti, Tf, Ts, fac,
								  0, 'o', 1, 'm',
								  randGen);
  }
  
  if (to_file == 1)	{
	outF = fopen(file_name_out, "w");
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
