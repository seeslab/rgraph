#include <stdio.h>
#include <math.h>
#include <stdlib.h>

long seed=42903;//6010;//42903;//20398;//82579;//9374;//

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"

#define rep 100

int main()
{
  // ----------------------------------------------------------------
  // Initialize the random number generator
  // ----------------------------------------------------------------
  gsl_rng *rand_gen;
  rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rand_gen, seed);

  // ----------------------------------------------------------------
  // Build the network---only use giant component
  // ----------------------------------------------------------------
  struct node_gra *total_net = NULL;
  struct node_gra *net = NULL;
  FILE *infile;

  infile = fopen("test.dat", "r");
  total_net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);
  net = GetLargestStronglyConnectedSet(total_net, 100000);

  // ----------------------------------------------------------------
  // Initialize the coclassification matrix
  // ----------------------------------------------------------------
  int i, j;
  int **coclas;
  int S = CountNodes(net);
  coclas = allocate_i_mat(S, S);
  for (i=0; i<S; i++)
    for (j=0; j<S; j++)
      coclas[i][j] = 0;

  // ----------------------------------------------------------------
  // Find the modules rep times and calculate coclassification matrix
  // ----------------------------------------------------------------
  struct group *part = NULL;
  struct node_gra *p1 = NULL;
  struct node_gra *p2 = NULL;
  for (i=0; i<rep; i++) {
    fprintf(stderr, "Repetition %d\n", i+1);
    part = SACommunityIdent(net,
			    0.0, -1.0, 0.0,
			    0.10,
			    S,         // #modules
			    'r',       // r=random initial conf
			    1,         // 0=No collective moves
			    'n',       // n=no output
			    rand_gen);
   
    // Update coclas matrix
    // --------------------------------------------------------------
    p1 = net;
    while ((p1 = p1->next) != NULL) {
      p2 = net;
      while ((p2 = p2->next) != NULL) {
	if (p1->inGroup == p2->inGroup)
	  coclas[p1->num][p2->num] += 1;
      }
    }

    // Remove the partition
    // --------------------------------------------------------------
    RemovePartition(part);
    part = NULL;
  }

  // ----------------------------------------------------------------
  // Output results
  // ----------------------------------------------------------------
  p1 = net;
  while ((p1 = p1->next) != NULL) {
    p2 = net;
    while ((p2 = p2->next) != NULL) {
      fprintf(stdout, "%s %s -1 -1 -1 -1 %lf\n",
	      p1->label, p2->label,
	      (double)coclas[p1->num][p2->num] / (double)rep);
    }
  }

  // ----------------------------------------------------------------
  // Free memory
  // ----------------------------------------------------------------
  RemoveGraph(total_net);
  RemoveGraph(net);
  free_i_mat(coclas,S);
  gsl_rng_free(rand_gen);

  return 0;
}
