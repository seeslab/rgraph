/*
  conflict.c
  $LastChangedDate: 2008-10-22 17:39:18 -0500 (Wed, 22 Oct 2008) $
  $Revision: 134 $
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "models.h"
#include "missing.h"
#include "recommend.h"
#include "conflict.h"

#define max(A, B) ((A > B)? A : B)

/*
  -----------------------------------------------------------------------------
  Return the score p(A_ij=1|A^O) of all pairs of nodes for a 2-state
  system. The ratings are a bipartite network with links
  (corresponding to observations) that have values 0 or 1.
  -----------------------------------------------------------------------------
*/
struct query **
AllLinkScore2State(struct binet *ratings,
			   int nIter,
			   gsl_rng *gen,
			   char verbose_sw,
			   int decorStep)
{
  int nnod1=CountNodes(ratings->net1), nnod2=CountNodes(ratings->net2);
  int nn1, nn2;
  struct node_gra *net1=NULL, *net2=NULL;
  struct group *part1=NULL, *part2=NULL;
  struct node_gra *p1=NULL, *p2=NULL, *node=NULL;
  struct node_gra **nlist1=NULL, **nlist2=NULL;
  struct group **glist1=NULL, **glist2=NULL;
  struct group *lastg=NULL;
  double H;
  int iter;
  int i, j;
  int **N1G2_0=NULL, **N2G1_0=NULL;
  int **N1G2_1=NULL, **N2G1_1=NULL;
  int **G1G2_0=NULL, **G2G1_0=NULL;
  int **G1G2_1=NULL, **G2G1_1=NULL;
  int LogChooseListSize = 500;
  double **LogChooseList=InitializeFastLogChoose(LogChooseListSize);
  int LogListSize = 5000;
  double *LogList=InitializeFastLog(LogListSize);
  int LogFactListSize = 2000;
  double *LogFactList=InitializeFastLogFact(LogFactListSize);
  struct node_lis *n1=NULL, *n2=NULL;
  double contrib;
  int norm = 0;
  int dice;
  int r, l;
  int ng1, ng2;
  int q, nq, nquery=nnod1*nnod2;
  struct query **querySet;

  fprintf(stderr, "Here\n");

  /*
    PRELIMINARIES
  */

  /* Create the queries and initialize scores */
  querySet = (struct query **) calloc(nquery, sizeof(struct query *));
  nq = 0;
  p1 = ratings->net1;
  while ((p1=p1->next) != NULL) {
    p2 = ratings->net2;
    while ((p2=p2->next) != NULL) {
      fprintf(stderr, "%s %s (%d/%d)\n", p1->label, p2->label, nq+1, nquery);
      querySet[nq] = CreateQuery(p1, p2);
      querySet[nq]->score = 0.0;
      nq++;
    }
  }

  fprintf(stderr, "Here\n");



  /* Map nodes and groups to a list for faster access */
  fprintf(stderr, ">> Mapping nodes and groups to lists...\n");
  nlist1 = (struct node_gra **) calloc(nnod1, sizeof(struct node_gra *));
  glist1 = (struct group **) calloc(nnod1, sizeof(struct group *));
  lastg = part1 = CreateHeaderGroup();
  p1 = net1 = ratings->net1;
  while ((p1 = p1->next) != NULL) {
    nlist1[p1->num] = p1;
    lastg = glist1[p1->num] = CreateGroup(lastg, p1->num);
  }
  nlist2 = (struct node_gra **) calloc(nnod2, sizeof(struct node_gra *));
  glist2 = (struct group **) calloc(nnod2, sizeof(struct group *));
  lastg = part2 = CreateHeaderGroup();
  p2 = net2 = ratings->net2;
  while ((p2 = p2->next) != NULL) {
    nlist2[p2->num] = p2;
    lastg = glist2[p2->num] = CreateGroup(lastg, p2->num);
  }

  /* Place each node in a group */
  fprintf(stderr, ">> Placing nodes in random partitions...\n");
  p1 = net1;
  ResetNetGroup(net1);
  while ((p1 = p1->next) != NULL) {
/*     dice = floor(gsl_rng_uniform(gen) * (double)nnod1); */
/*     AddNodeToGroup(glist1[dice], p1); */
    AddNodeToGroup(glist1[p1->num], p1);

  }
  p2 = net2;
  ResetNetGroup(net2);
  while ((p2 = p2->next) != NULL) {
/*     dice = floor(gsl_rng_uniform(gen) * (double)nnod2); */
/*     AddNodeToGroup(glist2[dice], p2); */
    AddNodeToGroup(glist2[p2->num], p2);
  }

  /* Get the initial group-to-group links matrix */
  fprintf(stderr, ">> Getting the initial group-to-group links matrix...\n");
  G1G2_0 = allocate_i_mat(nnod1, nnod2);
  G2G1_0 = allocate_i_mat(nnod2, nnod1);
  G1G2_1 = allocate_i_mat(nnod1, nnod2);
  G2G1_1 = allocate_i_mat(nnod2, nnod1);
  for (i=0; i<nnod1; i++) {
    for (j=0; j<nnod2; j++) {
      G1G2_0[i][j] = G2G1_0[j][i] =
	NWeightG2GLinks(glist1[i], glist2[j], (double)0);
      G1G2_1[i][j] = G2G1_1[j][i] =
	NWeightG2GLinks(glist1[i], glist2[j], (double)1);
    }
  }

  /* Get the initial node-to-group links matrix */
  fprintf(stderr, ">> Getting the initial node-to-group links matrix...\n");
  N1G2_0 = allocate_i_mat(nnod1, nnod2);
  N2G1_0 = allocate_i_mat(nnod2, nnod1);
  N1G2_1 = allocate_i_mat(nnod1, nnod2);
  N2G1_1 = allocate_i_mat(nnod2, nnod1);
  for (i=0; i<nnod1; i++) {
    for (j=0; j<nnod2; j++) {
      N1G2_0[i][j] = NWeightLinksToGroup(nlist1[i], glist2[j], (double)0);
      N2G1_0[j][i] = NWeightLinksToGroup(nlist2[j], glist1[i], (double)0);
      N1G2_1[i][j] = NWeightLinksToGroup(nlist1[i], glist2[j], (double)1);
      N2G1_1[j][i] = NWeightLinksToGroup(nlist2[j], glist1[i], (double)1);
    }
  }

  /* Get the initial number of non-empty groups */
  ng1 = NNonEmptyGroups(part1);
  ng2 = NNonEmptyGroups(part2);

  /*
    GET READY FOR THE SAMPLING
  */
  H = H2State(part1, part2);

  /* Get the decorrelation time */
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "# CALCULATING DECORRELATION TIME\n");
    fprintf(stderr, "# ------------------------------\n");
    break;
  }
  if (decorStep <= 0) {
    decorStep = GetDecorrelationStep2State(&H,
						   nlist1, nlist2,
						   glist1, glist2,
						   part1, part2,
						   nnod1, nnod2,
						   &ng1, &ng2,
						   N1G2_0, N2G1_0,
						   N1G2_1, N2G1_1,
						   G1G2_0, G2G1_0,
						   G1G2_1, G2G1_1,
						   LogList, LogListSize,
						   LogChooseList,
						   LogChooseListSize,
						   LogFactList, LogFactListSize,
						   gen, verbose_sw);
  }

  /* Thermalization */
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "#\n#\n# THERMALIZING\n");
    fprintf(stderr, "# ------------\n");
    break;
  }
  ThermalizeMC2State(decorStep, &H,
			     nlist1, nlist2, glist1, glist2,
			     part1, part2,
			     nnod1, nnod2, &ng1, &ng2,
			     N1G2_0, N2G1_0, N1G2_1, N2G1_1,
			     G1G2_0, G2G1_0, G1G2_1, G2G1_1,
			     LogList, LogListSize,
			     LogChooseList, LogChooseListSize,
			     LogFactList, LogFactListSize,
			     gen, verbose_sw);
  
  /*
    SAMPLIN' ALONG
  */
  switch (verbose_sw) {
  case 'd':
    break;
  default:
    H = 0; /* Reset the origin of energies to avoid huge exponentials */
    break;
  }
  for (iter=0; iter<nIter; iter++) {
    MCStep2State(decorStep, &H, nlist1, nlist2,
			 glist1, glist2, part1, part2,
			 nnod1, nnod2, &ng1, &ng2,
			 N1G2_0, N2G1_0, N1G2_1, N2G1_1,
			 G1G2_0, G2G1_0, G1G2_1, G2G1_1,
			 LogList, LogListSize,
			 LogChooseList, LogChooseListSize,
			 LogFactList, LogFactListSize,
			 gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf\n", iter, H, H2State(part1, part2));
      FPrintPartition(stderr, part1, 0);
      FPrintPartition(stderr, part2, 0);
      break;
    }

    /* Check if the energy has gone below a certain threshold and, if
       so, reset energies and start over */
    if (H < -400.0) {
      switch (verbose_sw) {
      case 'q':
	break;
      default:
	fprintf(stderr,
		"# System was not properly thermalized: starting over :(\n\n");
	break;
      }
      iter = 0;
      H = 0.0;
      norm = 0;
      for (q=0; q<nquery; q++) {
	querySet[q]->score = 0.0;
      }
    }

    /* Update normalization */
    norm += 1;

    /* Update the scores */
    for (q=0; q<nquery; q++) {
      l = G1G2_1[querySet[q]->n1->inGroup][querySet[q]->n2->inGroup];
      r = G1G2_0[querySet[q]->n1->inGroup][querySet[q]->n2->inGroup] +
	G1G2_1[querySet[q]->n1->inGroup][querySet[q]->n2->inGroup];
      contrib = (float)(l + 1) / (float)(r + 2);
      querySet[q]->score += contrib;
    }

  }  /* End of iter loop */

  /* Normalize the scores */
  for (q=0; q<nquery; q++)
    querySet[q]->score /= (double)norm;

  /* Free dynamically allocated memory */
  RemovePartition(part1);
  RemovePartition(part2);
  free(glist1);
  free(glist2);
  free(nlist1);
  free(nlist2);
  free_i_mat(G1G2_0, nnod1);
  free_i_mat(G2G1_0, nnod2);
  free_i_mat(G1G2_1, nnod1);
  free_i_mat(G2G1_1, nnod2);
  free_i_mat(N1G2_0, nnod1);
  free_i_mat(N2G1_0, nnod2);
  free_i_mat(N1G2_1, nnod1);
  free_i_mat(N2G1_1, nnod2);
  FreeFastLog(LogList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  FreeFastLogFact(LogFactList);

  /* Done */
  return querySet;
}
