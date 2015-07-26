#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <search.h>

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
  int i=0;
  int rolenb;
  fprintf (outf, "Label\tModule\tConnectivity\tParticipation\tRole\n");
  for (i=0;i<part->N;i++){
	rolenb = GetRole(participation[i],connectivity[i])+1;
	fprintf (outf, "%s\t%d\t%f\t%f\tR%d\n",
			 labels[i],
			 part->nodes[i]->module,
			 connectivity[i],participation[i],
			 rolenb);
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
AssignNodesToModulesFromFile(FILE *inF,
							 Partition *part,
							 char **labels){
  int mod;
  Node *node;
  char label[MAX_LABEL_LENGTH];
  char sep[2];
  int j = 0, nfields = 0, nnode = part->N;
  hcreate(part->N);
  int i;
  ENTRY e, *ep;
  
  for (i=0;i<nnode;i++){
	e.key = labels[i];
	e.data = (void *) i;
	ep = hsearch(e, ENTER);
  }

  while (!feof(inF)){
	nfields=fscanf(inF,"%[^\t\n]%[\t\n]",&label,&sep);
	if (nfields) {
	  e.key = label;
	  ep = hsearch(e, FIND);
	  i = ep->data;
	  nnode--;
	  
	  if (!part->modules[j]->size){
		part->nempty --;
		part->nodes[i]->module = j;
		part->modules[j]->size = 1;
		part->modules[j]->strength = part->nodes[i]->strength;
		part->modules[j]->first = part->nodes[i];
		part->modules[j]->last = part->nodes[i];
	  }
	  else{
		part->nodes[i]->module = j;
		part->modules[j]->size++;
		part->modules[j]->strength += part->nodes[i]->strength;
		part->modules[j]->last->next = part->nodes[i];
		part->nodes[i]->prev = part->modules[j]->last; 
		part->modules[j]->last = part->nodes[i];
	  }
	  if(sep[0]=='\n')
		j++;
	}
  }

  CompressPartition(part);
  hdestroy();
  return nnode;
}



int
main(int argc, char **argv)
{
  FILE *outF, *inF;
  int rgm;
  unsigned int seed = 1111;
  
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
  
  unsigned int N = 0, Ngroups=0, nochange_limit;
  
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
  int E;
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
  if (labels==NULL){
	perror("Error while setting labels");
	exit(1);
  }
  struct node_gra *node = net;
  while ((node=node->next)!=NULL)
	labels[node->num] = node->label;
  
  Ti = 1. / (double)N;
  Tf = 0.;
  nochange_limit = 25;
  fprintf(stderr, "# %d nodes to cluster.\n", N);
  if (!Ngroups)
	Ngroups = N;
  
  part = CreatePartition(N,Ngroups);

  if (!bipartite){
	// Allocate the memory.
	E = TotalNLinks(net, 1);
	adj = CreateAdjaArray(N,E);
	// Initialization.
	ComputeCost(net, adj, part);
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
	ComputeCostBipart(binet, adj, part, projected, weighted);
  }

  // SIMULATED ANNEALING CLUSTERING
  if (clustering){
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
  
  // Free memory
  if (bipartite)
	RemoveBipart(binet);
  else
  	RemoveGraph(net);
  if (roles){
	free(connectivity);
	free(participation);
  }
  free(labels);
  FreeAdjaArray(adj);
  FreePartition(part);  
  gsl_rng_free(randGen);
  return 0;
}
