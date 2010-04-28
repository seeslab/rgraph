/*
  main_graph1.c
  $LastChangedDate: 2009-02-02 18:27:16 -0600 (Mon, 02 Feb 2009) $
  $Revision: 158 $
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "graph.h"

#define EPS 1.e-6

int main()
{
  struct node_gra *net1=NULL, *net2=NULL, *netSum=NULL;
  FILE *inFile;
  double C, L;
  int wrong;

  /*
    ------------------------------------------------------------
    Build the networks
    ------------------------------------------------------------
  */
  inFile = fopen("test.dat", "r");
  net1 = FBuildNetwork(inFile, 0, 0, 0, 1);
  fclose(inFile);
  inFile = fopen("test_betw.dat", "r");
  net2 = FBuildNetwork(inFile, 0, 0, 0, 1);
  fclose(inFile);

  /*
    ------------------------------------------------------------
    Sum the networks
    ------------------------------------------------------------
  */
  netSum = AddTwoNetworks(net1, net2);

  /*
    ------------------------------------------------------------
    Verify the results and done
    ------------------------------------------------------------
  */
  C = ClusteringCoefficient(netSum);
  L = AverageInverseDistance(netSum);
  RemoveGraph(net1);
  RemoveGraph(net2);
  RemoveGraph(netSum);
  wrong = 0;
  if (fabs(C - 0.428571) > 1e-6) {
    wrong = 1;
    fprintf(stderr, "C: %lf != 0.428571 \tWRONG!\n", C);
  }
  else {
    fprintf(stderr, "C: %lf == 0.428571 \tOK!\n", C);
  }
  if (fabs(L - 0.448485) > 1e-6) {
    wrong = 1;
    fprintf(stderr, "L: %lf != 0.448485 \tWRONG!\n", L);
  }
  else {
    fprintf(stderr, "L: %lf == 0.448485 \tOK!\n", L);
  }
  if (wrong == 0)
    return 0;
  else
    return 1;
}
