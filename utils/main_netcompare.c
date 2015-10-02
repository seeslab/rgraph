/*
  main_netcompare.c
  $LastChangedDate$
  $Revision$
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"

int
main(int argc, char **argv)
{
  char *netFA, *netFB;
  FILE *inFile;
  struct node_gra *netA=NULL, *netB=NULL;
  int SA, SB, SAB;
  double p_nAB, p_nBA;
  int LA, LB, LAB;
  double p_lAB, p_lBA;

  /*
    ---------------------------------------------------------------------------
    Command line parameters
    ---------------------------------------------------------------------------
  */
  if (argc < 3) {
    printf("\nUse: netcompare.out net_file_A net_file_B\n\n");
    return -1;
  }
  netFA = argv[1];
  netFB = argv[2];

  /*
    ---------------------------------------------------------------------------
    Build the networks
    ---------------------------------------------------------------------------
  */
  inFile = fopen(netFA, "r");
  netA = FBuildNetwork(inFile, 0, 0, 0, 1);
  fclose(inFile);

  inFile = fopen(netFB, "r");
  netB = FBuildNetwork(inFile, 0, 0, 0, 1);
  fclose(inFile);

  /*
    ---------------------------------------------------------------------------
    Compare the networks
    ---------------------------------------------------------------------------
  */
  CompareTwoNetworks(netA, netB,
		     &SA, &SB, &SAB, &p_nAB, &p_nBA,
		     &LA, &LB, &LAB, &p_lAB, &p_lBA);
  
  printf("\nS_A = %d\n", SA);
  printf("S_B = %d\n", SB);
  printf("S_(A&B) = %d\n\n", SAB);

  printf("p_n(A|B) = %lf\n", p_nAB);
  printf("p_n(B|A) = %lf\n", p_nBA);
  
  printf("\nL_A = %d\n", LA);
  printf("L_B = %d\n", LB);
  printf("L_(A&B) = %d\n\n", LAB);

  printf("p_l(A|B) = %lf\n", p_lAB);
  printf("p_l(B|A) = %lf\n\n", p_lBA);


  /*
    ---------------------------------------------------------------------------
    Free memory and finish
    ---------------------------------------------------------------------------
  */
  RemoveGraph(netA);
  RemoveGraph(netB);
  return 0;
}
