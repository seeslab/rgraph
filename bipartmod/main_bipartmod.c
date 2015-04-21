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
  struct group *roles_part = NULL;
  struct node_gra *projected = NULL;
  gsl_rng *randGen;
  double degree_based = 0;
  double Ti, Tf;
  double Ts = 0.97;
  double fac = 1.0;
  double modularity = 0;
  char fn_array[256];
  char fno_array[256];
  char *file_name;
  char *file_name_out;
  int invert = 0;
  int weighted = 0;
  int output_type = 0;
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

	printf("\n# Use the strength (0) or degree (1) based role metrics: ");
	scanf("%d", &degree_based);


	printf("\n# Choose between tabular output (0), module partition (1) or roles partition (2):");
	scanf("%d", &output_type);}
  
  
  else{
	while ((c = getopt(argc, argv, "hwprmdf:s:i:c:o:")) != -1)
	  switch (c) {
	  case 'h':
		printf("\nUsage: bipartmod_cl [-f file] [-o file] [-s seed] [-i iter] [-c cool] [-pwmrh]\n"
			   "\nIf no arguments are provided, the program will fallback in interactive mode.\n\n"
			   "\t -f file: Name of the input network file (default: standard input)\n"
			   "\t -o file: Name of the output file (default: standard output)\n"
			   "\t -s seed: Random number seed (POSITIVE Integer, default 1111) \n"
			   "\t -i iter: Iteration factor (recommended 1.0, default 1.0)\n"
			   "\t -c cool: Cooling factor (recommended 0.950-0.995, default 0.97)\n "
			   "\t -p : Find modules for the second column (default: first) \n"
			   "\t -w : Read edge weights from the input's third column and uses the weighted modularity.\n"
			   "\t -d : Use degree based role metrics (default: strength based metrics)"
			   "\t -r : Output the roles partition rather than the default tabular output.\n"
			   "\t -m : Output the modules partition rather than the default tabular output. \n"
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
	  case 'm':
		output_type = 1;
		break;
	  case 'r':
		output_type = 2;
		break;
	  case 'o':
		to_file = 1;
		file_name_out = optarg;
		break;
	  case 'w':
		weighted = 1;
		break;
	  case 'd':
		degree_based = 1;
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

  // Compute the role partition if we have to.
  // Note that the roles are computed (but not as a partition)
  // if the tabular output is selected (output_type==0).
  if (output_type == 2){
	projected = ProjectBipart(binet);
	if (degree_based==1)
	  part = CatalogRoleIdent(projected,part);
	else
	  part = CatalogRoleIdentStrength(projected,part);
  }

  /*
    ------------------------------------------------------------
    Output
    ------------------------------------------------------------
  */

  if (to_file == 1)	{
	outF = fopen(file_name_out, "w");
	if (outF == NULL)
	  {
		printf("ERROR: Cannot write output (%s). \n", file_name_out);
		return(1);
	  }
	if (output_type != 0){
	  // Partition-type output (a la Netcarto).
	  FPrintPartition(outF, part, 0);
	  modularity = Modularity(part);
	  fprintf(outF, "# Modularity = %g\n", modularity);
	}
	else
	  // Tabular output. 
	  FPrintTabNodesBipart(outF, binet, part, degree_based);


	fclose(outF);
  }
  else{
	if (output_type != 0){
	  FPrintPartition(stdout, part, 0);
	  modularity = Modularity(part);
	  fprintf(stdout, "# Modularity = %g\n", modularity);
	}
	else
	  FPrintTabNodesBipart(stdout, binet, part, degree_based);
  }
  
  // Free memory
  RemovePartition(part);
  RemoveBipart(binet);
  if (output_type == 2)
	RemoveGraph(projected);
  gsl_rng_free(randGen);
}
