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
  char fn_array[256];
  char fno_array[256];
  char *file_name;
  char *file_name_out;
  int invert = 0;
  int weighted = 0;
  int roles = 0;
  extern char *optarg;
  extern int optind;
  int c; 
  
  /*
    ------------------------------------------------------------
    Arguments parsing
    ------------------------------------------------------------
  */
  
  if (argc == 1) {
	from_file = 1;
	to_file = 1;	
	file_name = fn_array;
	file_name_out = fno_array;
	
	printf("\n# Enter random number seed (POSITIVE integer): ");
	scanf("%d", &seed);
	printf("\n# Enter the name of the network file: ");
	scanf("%s", &fn_array);

	printf("\n# Enter the name of the output file: ");
	scanf("%s", &fno_array);

	printf("\n# Enter iteration factor (recommended 1.0): ");
	scanf("%lf", &fac);
  
	printf("\n# Enter the cooling factor (recommended 0.950-0.995): ");
	scanf("%lf", &Ts);

	printf("\n# Find modules from first column (0) or second column (1): ");
	scanf("%d", &invert);

	printf("\n# Use the weighted (0) or weighted (1) modularity. (Edge weight is extracted from the third column): ");
	scanf("%d", &weighted);

	printf("\n# Output the modules (0) or the roles (1) of the nodes. ");
	scanf("%d", &roles);
	
  }
  else{
	while ((c = getopt(argc, argv, "hwprf:s:i:c:o:")) != -1)
	  switch (c) {
	  case 'h':
		printf("\nUsage: bipartmod_cl [-f file] [-o file] [-s seed] [-i iter] [-c cool] [-pw]\n"
			   "\nIf no arguments are provided, the program will fallback in interactive mode.\n\n"
			   "\t -f file: Name of the input network file (default: standard input)\n"
			   "\t -o file: Name of the output file (default: standard output)\n"
			   "\t -s seed: Random number seed (POSITIVE Integer, default 1111) \n"
			   "\t -i iter: Iteration factor (recommended 1.0, default 1.0)\n"
			   "\t -c cool: Cooling factor (recommended 0.950-0.995, default 0.97)\n "
			   "\t -p : Find modules for the second column (default: first) \n"
			   "\t -w : Read edge weights from the input's third column and uses the weighted modularity.\n"
			   "\t -r : Output the roles rather than the modules.\n"
			   "\t -h : Display this message\n");
		return -1;
		break;
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
	  case 'r':
		roles = 1;
		break;
	  case 'o':
		to_file = 1;
		file_name_out = optarg;
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
	if (inF == NULL)
	  {
		printf("ERROR: No such file or directory (%s). \n", file_name);
		return(1);
	  }
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
	if (outF == NULL)
	  {
		printf("ERROR: Cannot write output (%s). \n", file_name_out);
		return(1);
	  }

    FPrintPartition(outF, part, 0);
	fclose(outF);
  }
  else{
	FPrintPartition(stdout, part, 0);
  }
  // Free memory
  // ------------------------------------------------------------
  RemovePartition(part);
  RemoveBipart(binet);
  gsl_rng_free(randGen);
}
