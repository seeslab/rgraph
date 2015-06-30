/*
  main_netcarto.c
  $LastChangedDate$
  $Revision$
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>
#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "sannealing.h"

int
main()
{
  long seed;
  FILE *inFile,*outFile;
  struct node_gra *net = NULL;
  struct node_gra *netran = NULL;
  struct node_gra *p = NULL;
  int S;
  struct group *part = NULL;
  struct group *roles = NULL;
  int i,t1;
  int rep;
  double realmod;
  double ranmodav, ranmodst;
  double *ranmodlis;
  double Tsched, Tf = 0.0;
  double iterfac;
  char netF[100];
  gsl_rng *rand_gen;


  /*
    ------------------------------------------------------------
    Prompt for user-defined parameters
    ------------------------------------------------------------
  */
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
