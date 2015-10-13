#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <search.h>

#include <gsl/gsl_rng.h>
#include "io.h"
#include "fillpartitions.h"
#include "sannealing.h"
#include "partition.h"

#define USAGE "Usage:\n\
\tnetcarto [-f FILE] [-o FILE] [-s SEED] [-i ITER] [-c COOL] [-wmr]\n\
\tnetcarto [-f FILE] [-o FILE] [-s SEED] [-i ITER] [-c COOL] [-wmr] -b [-t]\n\
\tnetcarto [-f FILE] [-o FILE] [-p FILE] [-w]\n\
\tnetcarto [-f FILE] [-o FILE] [-p FILE] [-w] -b [-t]\n\
\tnetcarto  -h\n"

#define ARGUMENTS "Arguments:\n\
\t -f FILE: Input network file name (default: '-', standard input),\n\
\t -o FILE: Output file name (default: '-', standard output),\n\
\t -s SEED: Random number generator seed (positive integer, default 1111),\n\
\t -i ITER: Iteration factor (recommended 1.0, default 1.0),\n\
\t -c COOL: Cooling factor (recommended 0.950-0.995, default 0.97),\n\
\t -p FILE: Partition file name to load and compute modularity and roles onto, \n\
\t -w : Read edge weights from the input's third column and uses the weighted modularity,\n\
\t -b : Use bipartite modularity,\n\
\t -r : Compute modularity roles,\n\
\t -t : [with -b only] Find modules for the second column (default: first),\n\
\t -h : Display this message.\n"

int
main(int argc, char **argv)
{
  FILE *outF, *inF;
  int rgm;
  unsigned int seed = 1111;
	unsigned int Ngroups = 0;
  struct binet *binet = NULL;
  struct node_gra *net = NULL;
  struct node_gra *projected = NULL;

  Partition *part = NULL;
  AdjaArray *adj = NULL;
  gsl_rng *randGen;

  double Ti, Tf;
  double Ts = 0.97;
  double fac = 1.0;
  double modularity = 0;
  double modularity_diag = 0;

  char fn_array[256];
  char fno_array[256];

  char *file_name;
  char *file_name_part;
  char *file_name_out;

  char **labels = NULL;

  unsigned int nochange_limit;
	unsigned int err;
  double proba_components = .5;

  int to_file = 0;
  int from_file = 0;
  int roles = 0;
  int invert = 0;
  int weighted = 0;
  int output_type = 0;
  int add_weight = 0;
  int louvain = 0;
  int bipartite = 0;
  int clustering = 1;
  int c;
  double *connectivity, *participation;
  int i;

  extern char *optarg;
  extern int optind;

  fprintf(stderr, "This is Netcarto.\n");

  //// Arguments parsing

  if (argc == 1) {
	printf(USAGE
		   "Use netcarto -h for more information.\n");
	return -1;
  }

  else{
	while ((c = getopt(argc, argv, "hbwtrmaf:s:i:c:o:S:p:C:")) != -1)
	  switch (c) {
	  case 'h':
		printf(USAGE ARGUMENTS);
		return -1;
		break;
	  case 'r':
		roles = 1;
		break;
	  case 'f':
		if(*optarg != '-')
		  from_file = 1;
		file_name = optarg;
		break;
	  case 'b':
		bipartite = 1;
		break;
	  case 'C':
		proba_components = atof(optarg);
		break;
	  case 's':
		seed = atoi(optarg);
		break;
	  case 'i':
		fac = atof(optarg);
		break;
	  case 'c':
		Ts = atof(optarg);
		break;
	  case 'p':
		file_name_part = optarg;
		clustering=0;
		break;
	  case 'm':
		output_type = 1;
		break;
	  case 'a':
		add_weight = 1;
		break;
	  case 'o':
		if(*optarg != '-')
		  to_file = 1;
		file_name_out = optarg;
		break;
	  case 'w':
		weighted = 1;
		break;
	  case 't':
		invert = 1;
		break;
	  }
  }

  //// Initialize the random number generator
  randGen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(randGen, seed);

  /*------------------------------------------------------------------
	Build the network
    ------------------------------------------------------------------ */
  fprintf(stderr, "Read input...\n");

  if (from_file == 1) {
  	inF = fopen(file_name, "r");
    if (inF == NULL){
    		printf("ERROR: No such file or directory (%s). \n", file_name);
    		return(1);
  	}
  }
  else	inF = stdin;

  unsigned int *nodes1 = NULL,*nodes2 = NULL;
  unsigned int E,N;
  double *weights = NULL;
	err = EdgeListFileInput(inF, weighted, bipartite+invert,
								 &nodes1, &nodes2, &weights, &labels,
								 &E, &N);
	if (Ngroups==0) Ngroups=N;

	if (!bipartite){
		part = CreatePartition(N,Ngroups);
		adj = CreateAdjaArray(N,E);
		err = EdgeListToAdjaArray(nodes1, nodes2,
															weights, adj, part, 1);
	}else{
		ProjectBipartEdgeList(nodes1, nodes2, weights, E,
		                      &part, &adj);
	}
  free(nodes1);
  free(nodes2);
  free(weights);
  if (from_file) fclose(inF);

  // SIMULATED ANNEALING CLUSTERING
  if (clustering){
    Ti = 1. / (double)N;
    Tf = 1e-200;
    nochange_limit = 25;
  	AssignNodesToModules(part,randGen);
  	GeneralSA(&part, adj, fac,
  			  Ti, Tf, Ts,
  			  proba_components, nochange_limit,
  			  randGen);
  	CompressPartition(part);
  }

  // READ PARTITION FROM A FILE.
  else {
	printf ("# Read partition from %s\n",file_name_part);
	inF = fopen(file_name_part, "r");
	if (inF == NULL){
		printf("ERROR: No such file or directory (%s). \n", file_name_part);
		return(1);
	}
	int i;
	i = AssignNodesToModulesFromFile(inF,part,labels);
	fclose(inF);
	if (i){
	  printf ("Error: %d nodes are missing from the partition file.\n",i);
	}
	printf ("# %d modules read \n",part->M);
  }

  modularity = PartitionModularity(part,adj,0);
  modularity_diag = PartitionModularity(part,adj,1);

  if(roles){
	connectivity = (double*) calloc(part->N,sizeof(double));
	participation = (double*) calloc(part->N,sizeof(double));
	if (connectivity==NULL || participation == NULL ){
	  perror("Error while setting roles metrics");
	  exit(1);
	}
	PartitionRolesMetrics(part, adj, connectivity, participation);
  }


  /* ------------------------------------------------------------
    Output
    ------------------------------------------------------------ */

  // Select output
  if (to_file == 1)	{
	outF = fopen(file_name_out, "w");
	if (outF == NULL) {
		printf("ERROR: Cannot write output (%s). \n", file_name_out);
		return(1);
	}
  } else outF = stdout;

  // Print output
  fprintf(outF,"# Modularity: %f\n",modularity);
  fprintf(outF,"# Modularity (with diagonal): %f\n",modularity_diag);

  if (roles)
	TabularOutput(outF, labels, part, connectivity, participation);
  else
	ClusteringOutput(outF, part, labels);

  // Close file if we need to.
  if (to_file == 1)	fclose(outF);

  //Free the memory
  if (roles){
	free(connectivity);
	free(participation);
  }
  for (i = 0; i < N; i++) {
    free(labels[i]);
  }
  free(labels);
  FreeAdjaArray(adj);
  FreePartition(part);
  gsl_rng_free(randGen);
  return 0;
}
