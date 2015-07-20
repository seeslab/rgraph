#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include "tools.h"
#include "louvain.c"
#include "graph.h"
#include "modules.h"
#include "bipartite.h"
#include "sannealing.h"
#include "partition.h"

#define USAGE "Usage:\n\
\tnetcarto [-f FILE] [-o FILE] [-s SEED] [-i ITER] [-c COOL] [-S NUM] [-wmr]\n\
\tnetcarto [-f FILE] [-o FILE] [-s SEED] [-i ITER] [-c COOL] [-S NUM] [-wmr] -b [-t]\n\
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
\t -S NUM: Number of graph shuffling to compute modularity's P-value (default 0),\n\
\t -T THRESHOLD: [with -L only] Louvain's method threshold (default is 0.00),\n\
\t -C PROBA: Proba of using connected components for the splits rather than SA. \n\
\t -w : Read edge weights from the input's third column and uses the weighted modularity,\n\
\t -b : Use bipartite modularity,\n\
\t -r : Compute modularity roles,\n\
\t -t : [with -b only] Find modules for the second column (default: first),\n\
\t -a : Sum the weights if an edges is defined several times (default: ignore subsequent definitions).\n\
\t -h : Display this message.\n"

void
TabularOutput(FILE *outf,
			  char **labels,
			  Partition *part,
			  double *connectivity,
			  double *participation){
  int i;
  printf ("Label\tModule\tConnectivity\tParticipation\tRole\n");
  for (i=0;i<part->N;i++){
	printf ("%s\t%d\t%f\t%f\tR%d\n",
			labels[i],
			part->nodes[i]->module,
			connectivity[i],participation[i],
			GetRole(participation[i],connectivity[i])+1);
  }
}

void
ClusteringOutput(FILE *outf,
				 Partition *part,
				 char **labels){
  int mod;
  Node *node;
  for(mod=0;mod<part->M;mod++){
	for(node = part->modules[mod]->first;node!=NULL; node=node->next){
	  fprintf(outf,"%s\t",labels[node->id]);
	}
	fprintf(outf,"\n");
  }
}


int
main(int argc, char **argv)
{
  FILE *outF, *inF;
  int rgm;
  unsigned int seed = 1111;
  int to_file = 0;
  int from_file = 0;
  struct binet *binet = NULL;
  struct node_gra *net = NULL;
  struct group *part = NULL;
  int roles = 0;
  struct node_gra *projected = NULL;
  gsl_rng *randGen;
  
  double degree_based = 0;
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
  unsigned int N = 0, Ngroups=0, nochange_limit;
  double proba_components = .5;
  int invert = 0;
  int weighted = 0;
  int output_type = 0;
  int add_weight = 0;
  int louvain = 0;
  int bipartite = 0;
  int clustering = 1;
  double louvain_threshold = 0;
  extern char *optarg;
  extern int optind;
  int c;
  fprintf(stderr, "This is Netcarto.\n");
  int E;
  double *connectivity, *participation;
  int i;

  //// Arguments parsing
  
  if (argc == 1) {
	printf(USAGE
		   "Use netcarto -h for more information.\n");
	return -1;
  }  
  
  else{
	while ((c = getopt(argc, argv, "hbwtrmdLaf:s:i:c:o:T:S:p:C:")) != -1)
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
	  case 'L':
		louvain = 1;
		break;
	  case 'T':
		louvain_threshold = atof(optarg);
		break;
	  case 'o':
		if(*optarg != '-')
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

  //// Initialize the random number generator
  randGen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(randGen, seed);

  /*------------------------------------------------------------------
	Build the network
    ------------------------------------------------------------------ */
  fprintf(stderr, "Building the network from input file...\n");  
  if (from_file == 1) {
	inF = fopen(file_name, "r");
	if (inF == NULL){
		printf("ERROR: No such file or directory (%s). \n", file_name);
		return(1);
	}
  }
  else
	inF = stdin;

  if (bipartite){
	binet = FBuildNetworkBipart(inF, weighted, add_weight); 
	if (invert == 1)
	  InvertBipart(binet);
	net = binet->net1;
  }
  else
	net = FBuildNetwork(inF,weighted,0,add_weight,1); // No autolink, symmetric.
  
  if (from_file)
	fclose(inF);
  
  if (binet == NULL && net == NULL){
	printf("Error reading input. \n");
	return(1);
  }
  
  /*
    ------------------------------------------------------------------
	Find the modules or load the partition.
    ------------------------------------------------------------------
  */
  N = CountNodes(net);
  labels = (char**) malloc(N*sizeof(char*));
  struct node_gra *node = net;
  while ((node=node->next)!=NULL)
	labels[node->num] = node->label;
  
  Ti = 1. / (double)N;
  Tf = 0.;
  nochange_limit = 25;
  fprintf(stderr, "# %d nodes to cluster.\n", N);
  if (!Ngroups)
	Ngroups = N;
  Partition *p = NULL;
  AdjaArray *adj = NULL;
  p = CreatePartition(N,Ngroups);

  if (!bipartite){
	// Allocate the memory.
	E = TotalNLinks(net, 1);
	adj = CreateAdjaArray(N,E);
	// Initialization.
	ComputeCost(net, adj, p);
  }
  else{   
	if (!weighted)
	  projected = ProjectBipart(binet);
	else
	  projected = ProjectBipartWeighted(binet);

	// Allocate the memory for the static graph.
	E = TotalNLinks(projected, 1);
	adj = CreateAdjaArray(N,E);
	// Initialization of the static graph.
	ComputeCostBipart(binet, adj, p, projected, weighted);
  }
  if (clustering){
	AssignNodesToModules(p,randGen);
	p = GeneralSA(p, adj, fac,
				  Ti, Tf, Ts,
				  proba_components, nochange_limit,
				  randGen);
	CompressPartition(p);
  } else{
	printf ("# Read partition from %s\n",file_name_part);
	inF = fopen(file_name_part, "r");
	if (inF == NULL){
		printf("ERROR: No such file or directory (%s). \n", file_name_part);
		return(1);
	}
	fclose(inF);
	printf ("# Not implemented yet ! \n");
  }
  
  modularity = PartitionModularity(p,adj,0);
  modularity_diag = PartitionModularity(p,adj,1);
  
  if(roles){
	connectivity = (double*) malloc(p->N*sizeof(double));
	participation = (double*) malloc(p->N*sizeof(double));
	PartitionRolesMetrics(p,adj, connectivity, participation);
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
	TabularOutput(outF, labels, p, connectivity, participation);
  else
	ClusteringOutput(outF, p, labels);

  // Close file if we need to. 
  if (to_file == 1)	fclose(outF);
  
  // Free memory
  free(labels);
  if (bipartite)
	RemoveBipart(binet);
  else
	RemoveGraph(net);
  gsl_rng_free(randGen);
}
