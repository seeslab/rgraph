/*
  main_bipart2.c
  $LastChangedDate: 2007-06-29 14:58:36 -0500 (Fri, 29 Jun 2007) $
  $Revision: 70 $
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "graph.h"
#include "bipartite.h"

int main()
{
  FILE *inf = NULL;
  struct binet *binet = NULL;
  struct group *part = NULL;
  void *labelDict = NULL;
  struct node_gra *target = NULL;
  double P;
  int result = 0;

  /*
    ------------------------------------------------------------
    Build the network and the partition
    ------------------------------------------------------------
  */
  fprintf(stderr, "Creating the network...\n");
  inf = fopen("testbinet.dat", "r");
  binet = FBuildNetworkBipart(inf, 0, 0);
  fclose(inf);

  fprintf(stderr, "Creating the partition...\n");
  inf = fopen("testbinetpart.dat", "r");
  part = FCreatePartition(inf);
  fclose(inf);

  fprintf(stderr, "Mapping the partition...\n");
  MapPartToNet(part, binet->net1);

  /*
    ------------------------------------------------------------
    Network properties
    ------------------------------------------------------------
  */
  fprintf(stderr, "Calculating participation...\n");
  labelDict = MakeLabelDict(binet->net1);
  target = GetNodeDict("1", labelDict);
  P = ParticipationCoefficientBipart(target);
  fprintf(stderr, "P(%s) = %lf\n", target->label, P);
  if (fabs(P - 12./25.) > 1.e-10)
    result = 1;  /* Wrong result */

  target = GetNodeDict("5", labelDict);
  P = ParticipationCoefficientBipart(target);
  fprintf(stderr, "P(%s) = %lf\n", target->label, P);
  if (fabs(P - 1./2.) > 1.e-10)
    result = 1;  /* Wrong result */

  target = GetNodeDict("2", labelDict);
  P = ParticipationCoefficientBipart(target);
  fprintf(stderr, "P(%s) = %lf\n", target->label, P);
  if (fabs(P) > 1.e-10)
    result = 1;  /* Wrong result */

  /* Participation coefficient for bipartite networks */
  double m, s, mi, ma;
  StatisticsParticipationCoefficientBipart(binet->net1,
					   &m,
					   &s,
					   &mi,
					   &ma);
  if (fabs(m - 1.48/6.) > 1.e-10)
    result = 1;  /* Wrong result */
  if (fabs(s - 0.24675674031) > 1.e-10)
    result = 1;  /* Wrong result */
  if (fabs(mi - 0.0) > 1.e-10)
    result = 1;  /* Wrong result */
  if (fabs(ma - 0.5) > 1.e-10)
    result = 1;  /* Wrong result */

  /* Group size statistical properties */
  GroupSizeStatistics(part, &m, &s, &mi, &ma);
  if (fabs(m - 3.0) > 1.e-10)
    result = 1;  /* Wrong result */
  if (fabs(s - 1.0) > 1.e-10)
    result = 1;  /* Wrong result */
  if (fabs(mi - 2.0) > 1.e-10)
    result = 1;  /* Wrong result */
  if (fabs(ma - 4.0) > 1.e-10)
    result = 1;  /* Wrong result */


  /*
    ------------------------------------------------------------
    Free memory
    ------------------------------------------------------------
  */
  RemoveBipart(binet);
  FreeLabelDict(labelDict);

  return result;
}
