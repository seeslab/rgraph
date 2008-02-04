/*
  main_netcarto_cl.c
  $LastChangedDate: 2007-09-13 10:24:20 -0500 (Thu, 13 Sep 2007) $
  $Revision: 94 $

  Similar to netcarto, but arguments are given in the command line
  (plus initial temperature can be freely specified).
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prng.h"
#include "tools.h"
#include "graph.h"
#include "modules.h"

int
main(int argc, char **argv)
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
  double Tsched, Ti, Tf = -1.0;
  double iterfac;
  char *netF;
  struct prng *rand_gen;


  /*
    ------------------------------------------------------------
    Prompt for user-defined parameters
    ------------------------------------------------------------
  */
  if (argc < 7) {
    printf("\nUse: netcarto_cl net_file_name seed T_ini iteration_factor");
    printf(" cooling_factor randomizations\n\n");
    return -1;
  }

  netF = argv[1];
  seed = atoi(argv[2]);
  Ti = atof(argv[3]);
  iterfac = atof(argv[4]);
  Tsched = atof(argv[5]);
  rep = atoi(argv[6]);
  
  /* Default values */
  if (iterfac < 0.0)
    iterfac = 1.0;
  if (Tsched < 0.0)
    Tsched = 0.995;

  /*
    ------------------------------------------------------------
    Initialize the random number generator and some vectors
    ------------------------------------------------------------
  */
  rand_gen = prng_new("mt19937(1111)");
  prng_seed(rand_gen, seed);
  ranmodlis = allocate_d_vec(rep);

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
  if (Ti < 0.0)
    Ti = 2.0 / (double)S; /* Default initial temperature */
  fprintf(stderr, "\n# The network has %d nodes\n", S);

  /*
    ------------------------------------------------------------
    Find and print the modules
    ------------------------------------------------------------
  */
  fprintf(stderr, "\n# Starting the simulated annealing\n");
  fprintf(stderr, "\n# 1/T\tM\tT\n");

  part = SACommunityIdent(net,
			  Ti, Tf, Tsched,
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
	    p->label, CountLinks(p),
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
  free_d_vec(ranmodlis);
  RemovePartition(part);
  RemovePartition(roles);
  RemoveGraph(net);
  if (netran != NULL)
    RemoveGraph(netran);
  
  return 0;
}
