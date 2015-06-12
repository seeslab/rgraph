/*
  main_countlinks.c
  $LastChangedDate: 2008-10-07 16:15:28 -0500 (Tue, 07 Oct 2008) $
  $Revision: 129 $
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "graph.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL;
  struct node_gra *net=NULL;
  int N, L;

  /*
    ---------------------------------------------------------------------------
    Command line parameters
    ---------------------------------------------------------------------------
  */
  if (argc < 2) {
    printf("\nUse: countlinks.out net_file\n\n");
    return -1;
  }
  netF = argv[1];

  /*
    ---------------------------------------------------------------------------
    Build the network
    ---------------------------------------------------------------------------
  */
  infile = fopen(netF, "r");
  net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);

  /*
    ---------------------------------------------------------------------------
    Output results and finish
    ---------------------------------------------------------------------------
  */
  N = CountNodes(net);
  L = TotalNLinks(net, 1);
  printf("Nodes: %d\tLinks: %d\t<k> = %g\n",
	 N, L, (double)(2 * L) / (double)N);

  RemoveGraph(net);
  return 0;
}
