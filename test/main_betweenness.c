#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "graph.h"

#define EPS 1.e-6

int
main()
{
  FILE *infile=NULL;
  struct node_gra *net=NULL;
  struct node_gra *p=NULL;
  double rightResult[8];

  /* The expected results */
  rightResult[1] = 4 + 2. / 3.;
  rightResult[2] = 8 + 1. / 3.;
  rightResult[3] = 3 + 1. / 3.;
  rightResult[4] = 2 + 1. / 3.;
  rightResult[5] = 2 + 1. / 3.;
  rightResult[6] = 4 + 2. / 3.;
  rightResult[7] = 8 + 1. / 3.;

  /* Build the network */
  infile = fopen("test_betw.dat", "r");
  net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);

  /* Calculate betweenness */
  CalculateNodeBetweenness(net);

  /* Some statistics */
  double betMean, betStddev, betMin, betMax;
  NodeBetweennessStatistics(net, &betMean, &betStddev, &betMin, &betMax);
  printf("b_mean\t= %g\n", betMean);
  printf("b_stdd\t= %g\n", betStddev);
  printf("b_min\t= %g\n", betMin);
  printf("b_max\t= %g\n", betMax);

  /* Check results */
  p = net;
  while ((p = p->next) != NULL) {
    printf("%s\t%g %g ", p->label, p->dvar1, rightResult[atoi(p->label)]);
    if (p->dvar1 - rightResult[atoi(p->label)] > EPS) {
      printf("WRONG!\n");
      return 1;
    }
    else
      printf("OK\n");
  }

  /* Free memory */
  RemoveGraph(net);

  /* Done */
  return 0;
}
