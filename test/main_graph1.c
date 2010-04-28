/*
  main_graph1.c
  $LastChangedDate$
  $Revision$
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include "graph.h"
#include "models.h"
#include "modules.h"

#define EPS 1.e-6

int main()
{
  FILE *infile=NULL;
  struct node_gra *net = NULL;
  int repeat = 1;
  gsl_rng *randGen;

  /*
    ------------------------------------------------------------
    Build the network
    ------------------------------------------------------------
  */
  infile = fopen("test_graph1.dat", "r");
  net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);
  FPrintPajekFile("graph1.net", net, 0, 0, 1);

  /*
    --------------------------------------------------------------
    Initialize the random number generator
    --------------------------------------------------------------
  */
  randGen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(randGen, 3333);

  /*
    ------------------------------------------------------------
    Network properties
    ------------------------------------------------------------
  */
  double dtemp;
  int itemp;
  printf("S = %d\n", CountNodes(net));

  dtemp = ClusteringCoefficient(net);
  if (fabs(dtemp - 0.1875) > EPS)
    return 1;
  else
    printf("C = %g\tOK\n", dtemp);

  printf("A = %g\n", Assortativity(net));

  dtemp = AverageInverseDistance(net);
  if (fabs(dtemp - 19. / 45.) > EPS) {
    printf("P = %g\tWRONG!!\n", dtemp);
    return 1;
  }
  else
    printf("P = %g\tOK\n", dtemp);

  struct node_gra *n1 = net->next;
  struct node_gra *n2 = n1->next;
  struct node_gra *n4 = n2->next->next;
  struct node_gra *n5 = n4->next;
  struct node_gra *n8 = n5->next->next->next;
  dtemp = JaccardIndex(n1, n2);
  if (fabs(dtemp - 0.25) > EPS)
    return 1;
  else
    printf("J[1][2] = %g\tOK\n", dtemp);
  dtemp = JaccardIndex(n2, n4);
  if (fabs(dtemp - 0.50) > EPS)
    return 1;
  else
    printf("J[2][4] = %g\tOK\n", dtemp);
  dtemp = JaccardIndex(n5, n8);
  if (fabs(dtemp - 0.50) > EPS)
    return 1;
  else
    printf("J[5][8] = %g\tOK\n", dtemp);
  
  if ((itemp = CommonNeighbors(n1, n2)) != 1)
    return 1;
  else
    printf("C[1][2] = %d\tOK\n", itemp);
  if ((itemp = CommonNeighbors(n2, n4)) != 1)
    return 1;
  else
    printf("C[2][4] = %d\tOK\n", itemp);
  if ((itemp = CommonNeighbors(n5, n8)) != 2)
    return 1;
  else
    printf("C[5][8] = %d\tOK\n", itemp);
      
  /*
    ------------------------------------------------------------
    Partition
    ------------------------------------------------------------
  */
  struct group *part = NULL;
  part = SACommunityIdent(net, 2. / CountNodes(net), 0.0, 0.995,
			  2.0, 2, 'o', 1, 'n', randGen);
  MapPartToNet(part, net);
  printf("M = %g\n", Modularity(part));
  RemovePartition(part);

  /*
    ------------------------------------------------------------
    Free memory
    ------------------------------------------------------------
  */
  RemoveGraph(net);
  gsl_rng_free(randGen);

  return 0;
}
