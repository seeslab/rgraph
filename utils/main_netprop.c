#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* #include "prng.h" */
#include "tools.h"
#include "graph.h"
#include "modules.h"
/* #include "missing.h" */
#include "matrix.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL;
  struct node_gra *net=NULL;
  struct node_gra *components[1000];
  int i,ncomponents,sizeGC=0,sizeothers=0;
/*   struct prng *rand_gen; */
/*   int seed; */

  /* Command line parameters */
  if (argc < 2) {
    printf("\nUse: netprop net_file\n\n");
    return -1;
  }
  netF = argv[1];
/*   seed = atoi(argv[2]); */
/*   rand_gen = prng_new("mt19937(1111)"); */
/*   prng_seed(rand_gen, seed); */

  /* Build the network */
  infile = fopen(netF, "r");
  net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);

  /*
    Calculate and print a variety of network properties
  */
  /* Calculate the size of the giant component and the number of components */
  ncomponents = GetAllConnectedSets(net,components);
  for(i=0;i<ncomponents;i++)
  {
    if (CountNodes(components[i])>sizeGC) sizeGC = CountNodes(components[i]);
    sizeothers += CountNodes(components[i]);
  }

  printf("Num_components %d\n", ncomponents);
  printf("Size_GC %d\n", sizeGC);
  if (ncomponents>1) printf("Av_size_other %g\n", (sizeothers-sizeGC)/(ncomponents-1.0));
  else printf("Av_size_other 0.0\n");

  /* Number of nodes, links, and so on */
  printf("Nodes %d\n", CountNodes(net));
  printf("Links %d\n", TotalNLinks(net, 1));
  printf("Av_degree %g\n", AverageDegree(net, 1));

  /* Distance, clustering... */
  printf("Inverse_path_length %lf\n", AverageInverseDistance(net));
  printf("Clustering %lf\n", ClusteringCoefficient(net));
  printf("Assortativity %lf\n", Assortativity(net));

  /* Betweenness */
  double betMean, betStddev, betMin, betMax;
  NodeBetweennessStatistics(net, &betMean, &betStddev, &betMin, &betMax);
  printf("Betweenness_mean %lf\n", betMean);
  printf("Betweenness_std %lf\n", betStddev);
  printf("Betweenness_min %lf\n", betMin);
  printf("Betweenness_max %lf\n", betMax);

  /* Synchronizability */
  printf("Synchronizability %lf\n", Synchronizability(net));

  printf("K_K2 %.10lf\n",  AverageDegree(net, 1)/ AverageSquaredDegree(net));

/*   /\* Modularity *\/ */
/*   struct group *part; */
/*   part = SACommunityIdent(net, */
/* 			  2.0 / (double)CountNodes(net), 0.0, */
/* 			  0.995, 1.0, 0, 'o', 1, 'q', rand_gen); */
/*   printf("Modularity %lf\n", Modularity(part)); */
/*   RemovePartition(part); */

  /* Finish */
  RemoveGraph(net);
  return 0;
}
