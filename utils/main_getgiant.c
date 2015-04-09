#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "tools.h"
#include "graph.h"

int
main(int argc, char **argv)
{
  char *netF;
  FILE *infile=NULL;
  struct node_gra *net=NULL, *giant=NULL;

  /* Command line parameters */
  if (argc < 2) {
    printf("\nUse: getgiant net_file\n\n");
    return -1;
  }
  netF = argv[1];

  /* Build the network */
  infile = fopen(netF, "r");
  net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);

  /* Get the giant component and print it */
  giant = GetLargestWeaklyConnectedSet(net, 10000000);
  FPrintNetAdjacencyList(stdout, giant, 0, 1);

  /* Finish */
  RemoveGraph(net);
  RemoveGraph(giant);
  return 0;
}
