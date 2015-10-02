#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
 
#include <gsl/gsl_rng.h>
#include "tools.h"
#include "graph.h"
#include "modules.h"

#define USAGE "Usage:\n\
\tnetcarto_legacy [-f FILE] [-s SEED] [-i ITER] [-c COOL] [-S SHUF]\n\
\tnetcarto_legacy  -h\n"

#define ARGUMENTS "Arguments:\n\
\t -f FILE: Input network file name (default: '-', standard input),\n\
\t -s SEED: Random number generator seed (positive integer, default 1111),\n\
\t -i ITER: Iteration factor (recommended 1.0, default 1.0),\n\
\t -c COOL: Cooling factor (recommended 0.950-0.995, default 0.97),\n\
\t -r NUM: Number of graph shuffling to compute modularity's P-value (default 0),\n\
\t -h : Display this message.\n"

int
main(int argc, char **argv)
{
  long seed = 1111;
  FILE *inFile,*outFile;
  struct node_gra *net = NULL;
  struct node_gra *netran = NULL;
  struct node_gra *p = NULL;
  int S;
  int c;
  struct group *part = NULL;
  struct group *roles = NULL;
  int i,t1;
  int rep = 0;
  double realmod;
  double ranmodav, ranmodst;
  double *ranmodlis;
  double Tsched = 0.995, Tf = 0.0;
  double iterfac = 1.0;
  gsl_rng *rand_gen;
  char *netF;
  /*
    ------------------------------------------------------------
    Prompt for user-defined parameters
    ------------------------------------------------------------
  */
  printf ("%d\n",argc);
  if (argc == 1) {

  printf("\n# Enter random number seed (POSITIVE integer): ");
  scanf("%d", &seed);

  printf("\n# Enter the name of the network file: ");
  scanf("%s", &netF);

  printf("\n# Enter iteration factor (recommended 1.0): ");
  scanf("%lf", &iterfac);

  printf("\n# Enter the cooling factor (recommended 0.950-0.995): ");
  scanf("%lf", &Tsched);

  printf("\n# Enter the number of randomizations: ");
  scanf("%d", &rep);
  } else {
while ((c = getopt(argc, argv, "hf:s:i:c:r:")) != -1)
	  switch (c) {
	  case 'h':
		printf(USAGE ARGUMENTS);
		return -1;
		break;
	  case 'f':
		netF = optarg;
		break;
	  case 's':
		seed = atoi(optarg);
		break;
	  case 'r':
		rep = atoi(optarg);
		break;
	  case 'i':
		iterfac = atof(optarg);
		break;
	  case 'c':
		Tsched = atof(optarg);
		break;
	  }
  }
  ranmodlis = allocate_d_vec(rep);

  /*
    ------------------------------------------------------------
    Initialize the random number generator
    ------------------------------------------------------------
  */
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, seed);

  /*
    ------------------------------------------------------------
    Build the network
    ------------------------------------------------------------
  */
  fprintf(stderr, "\n# Creating the network\n");

  inFile=fopen(netF,"r");
  net = FBuildNetwork(inFile, 0, 0, 0, 1);
  fclose(inFile);

  S = CountNodes(net);
  fprintf(stderr, "\n# The network has %d nodes\n", S);

  /*
    ------------------------------------------------------------
    Find and print the modules
    ------------------------------------------------------------
  */
  fprintf(stderr, "\n# Starting the simulated annealing\n");
  fprintf(stderr, "\n# 1/T\tM\tT\n");

  part = SACommunityIdent(net,
			  2.0 / (double)S, Tf, Tsched,
			  iterfac, 0, 'o', 1, 's', rand_gen);

  outFile = fopen("modules.dat", "w");
  FPrintPartition(outFile, part, 0);
  realmod = Modularity(part);
  fprintf(outFile, "# Modularity = %g\n", realmod);
  fclose(outFile);

  FPrintPajekFile("network.net", net, 0, 0, 1);
  FPrintPajekPartitionFile("modules.clu", net);

  /*
    ------------------------------------------------------------
    Calculate and print node properties
    ------------------------------------------------------------
  */
  outFile=fopen("node_prop.dat","w");
  p = net;
  while ((p = p->next) != NULL) {
    fprintf(outFile,"%s %d %lf %lf\n",
	    p->label, NodeDegree(p),
	    ParticipationCoefficient(p),
	    WithinModuleRelativeDegree(p, part));
  }
  fclose(outFile);

  /*
    ------------------------------------------------------------
    Find and print the roles
    ------------------------------------------------------------
  */
  fprintf(stderr, "\n# Finding node roles\n");
  roles = CatalogRoleIdent(net, part);

  outFile = fopen("roles.dat","w");
  FPrintPartition(outFile, roles, 0);
  fclose(outFile);

  MapPartToNet(roles, net);
  FPrintPajekPartitionFile("roles.clu", net);

  /*
    ------------------------------------------------------------
    Calculate properties for the randomized graph
    ------------------------------------------------------------
  */
  ranmodav = 0.0;
  for (i=0; i<rep; i++) {
    fprintf(stderr, "\nRandomization %d\n", i+1);

    RemovePartition(part);
    if ( netran != NULL)
      RemoveGraph(netran);

    netran = RandomizeSymmetricNetwork(CopyNetwork(net), 100, rand_gen);
    part = SACommunityIdent(netran,
			    2.0 / (double)S, Tf, Tsched,
			    iterfac, 0, 'o', 1, 'n', rand_gen);
    ranmodlis[i] = Modularity(part);
    ranmodav += ranmodlis[i];
  }

  /* Average and standard deviation */
  ranmodav /= (double)rep;
  ranmodst = 0.0;
  for (i=0; i<rep; i++) {
    ranmodst += (ranmodlis[i] - ranmodav) * (ranmodlis[i] - ranmodav);
  }
  ranmodst = sqrt(ranmodst / (double)rep);

  outFile = fopen("randomized_mod.dat", "w");
  fprintf(outFile, "# M\tMrand\tsimga_Mrand\n");
  fprintf(outFile, "%lf %lf %lf\n", realmod, ranmodav, ranmodst);
  fclose(outFile);

  /*
    ------------------------------------------------------------
    Free memory
    ------------------------------------------------------------
  */
  gsl_rng_free(rand_gen);
  free_d_vec(ranmodlis);
  RemovePartition(part);
  RemovePartition(roles);
  RemoveGraph(net);
  if (netran != NULL)
    RemoveGraph(netran);
  return 0;
}
