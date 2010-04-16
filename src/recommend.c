/*
  recommend.c
  $LastChangedDate: 2008-10-22 17:39:18 -0500 (Wed, 22 Oct 2008) $
  $Revision: 134 $
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "prng.h"

#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "models.h"
#include "missing.h"
#include "recommend.h"

#define max(A, B) ((A > B)? A : B)

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Auxiliary functions
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
struct query *
CreateQuery(struct node_gra *node1, struct node_gra *node2)
{
  struct query *theQuery; 
  theQuery = (struct query *)calloc(1, sizeof(struct query));
  theQuery->n1 = node1;
  theQuery->n2 = node2;
  return theQuery;
}

void
FreeQuery(struct query *q)
{
  free(q);
  return;
}

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Recommender functions
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  -----------------------------------------------------------------------------
  Hamiltoninan H for 2-state recommender. query_list contains the list
  of unobserved links
  -----------------------------------------------------------------------------
*/
double
H2State(struct group *part1, struct group *part2,
	struct query **query_list, int nquery)
{
  struct group *g1=NULL, *g2=NULL;
  double r, l, H=0.0;
  int q;

  g1 = part1;
  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
      g2 = part2;
      while ((g2 = g2->next) != NULL) {
	if (g2->size > 0) {
	  r = g1->size * g2->size;
	  l = NG2GLinks(g1, g2);
	  /* Discount unobserved query links */
	  for (q=0; q<nquery; q++) {
	    if ((query_list[q]->n1)->inGroup == g1->label &&
		(query_list[q]->n2)->inGroup == g2->label) {
	      r--;
	      l -= IsThereLink(query_list[q]->n1, query_list[q]->n2);
	    }
	  }
/* 	  fprintf(stderr, "%d-%d: %g/%g\n", g1->label+1, g2->label+1, l, r); */
	  H += log(r + 1) + LogChoose(r, l);
	}
      }
    }
  }

  return H;
}

/*
  -----------------------------------------------------------------------------
  Do a Monte Carlo step for the 2-state recommender system.
  -----------------------------------------------------------------------------
*/
void
MCStep2State(int factor,
	     double *H,
	     struct query **query_list, int nquery, 
	     struct node_gra **nlist1, struct node_gra **nlist2,
	     struct group **glist1, struct group **glist2,
	     struct group *part1, struct group *part2,
	     int nnod1, int nnod2,
	     int **G1G2, int **G2G1,
	     int *n2gList,
	     double **LogChooseList,
	     int LogChooseListSize,
	     struct prng *gen)
{
  double dH;
  struct group *oldg, *newg, *g, *g2;
  int ***G2G, ***G2Ginv;
  int move;
  double set_ratio;
  struct node_gra *node=NULL;
  int dice, r, l;
  int oldgnum, newgnum, q;
  int i, nnod, ngroup;
  int j; 
  int move_in;

  /* Ratio of moves in each of the sets */
  set_ratio = (double)(nnod1*nnod1-1) / (double)(nnod1*nnod1+nnod2*nnod2-2);

  for (move=0; move<(nnod1+nnod2)*factor; move++) {
    /* The move */
    if (prng_get_next(gen) < set_ratio) { /* move in first set */
      move_in = 1;
      G2G = &G1G2;
      G2Ginv = &G2G1;
      nnod = nnod1;
      ngroup = nnod2;
      dice = floor(prng_get_next(gen) * (double)nnod1);
      node = nlist1[dice];
      oldgnum = node->inGroup;
      do {
	newgnum = floor(prng_get_next(gen) * (double)nnod1);
      } while (newgnum == oldgnum);
      oldg = glist1[oldgnum];
      newg = glist1[newgnum];
      g = g2 = part2;
    }
    else {                                /* move in second set */
      move_in = 2;
      G2G = &G2G1;
      G2Ginv = &G1G2;
      nnod = nnod2;
      ngroup = nnod1;
      dice = floor(prng_get_next(gen) * (double)nnod2);
      node = nlist2[dice];
      oldgnum = node->inGroup;
      do {
	newgnum = floor(prng_get_next(gen) * (double)nnod2);
      } while (newgnum == oldgnum);
      oldg = glist2[oldgnum];
      newg = glist2[newgnum];
      g = g2 = part1;
    }

    /* The change of energy */
    /* Old configuration contribution */
    dH = 0.0;
    while ((g=g->next) != NULL) {
      if (g->size > 0) {  /* group is not empty */
      	n2gList[g->label] = NLinksToGroup(node, g);
	/* old configuration, old group */
	r = oldg->size * g->size;
	l = (*G2G)[oldgnum][g->label];
	for (q=0; q<nquery; q++) {
	  if ((move_in == 1 &&
	       query_list[q]->n1->inGroup == oldg->label &&
	       query_list[q]->n2->inGroup == g->label) ||
	      (move_in == 2 &&
	       query_list[q]->n1->inGroup == g->label &&
	       query_list[q]->n2->inGroup == oldg->label)) {
	    r--;
	    l -= IsThereLink(query_list[q]->n1, query_list[q]->n2);
	  }
	}
	dH -= log(r + 1) + LogChoose(r, l);
	/* old configuration, new group */
	r = newg->size * g->size;
	l = (*G2G)[newgnum][g->label];
	for (q=0; q<nquery; q++) {
	  if ((move_in == 1 &&
	       query_list[q]->n1->inGroup == newg->label &&
	       query_list[q]->n2->inGroup == g->label) ||
	      (move_in == 2 &&
	       query_list[q]->n1->inGroup == g->label &&
	       query_list[q]->n2->inGroup == newg->label)) {
	    r--;
	    l -= IsThereLink(query_list[q]->n1, query_list[q]->n2);
	  }
	}
	dH -= log(r + 1) + LogChoose(r, l);
      }
      else { /* group is empty */
	n2gList[g->label] = 0;
      }
    }

    /* Tentatively move the node to the new group and update G2G matrices */
    MoveNode(node, oldg, newg);
    for (i=0; i<ngroup; i++) {
      (*G2G)[oldgnum][i] -= n2gList[i];
      (*G2G)[newgnum][i] += n2gList[i];
      (*G2Ginv)[i][oldgnum] -= n2gList[i];
      (*G2Ginv)[i][newgnum] += n2gList[i];
    }

    /* New configuration contribution */
    while ((g2=g2->next) != NULL) {
      if (g2->size > 0) {  /* group is not empty */
	/* new configuration, old group */
	r = oldg->size * g2->size;
	l = (*G2G)[oldgnum][g2->label];
	for (q=0; q<nquery; q++) {
	  if ((move_in == 1 &&
	       query_list[q]->n1->inGroup == oldg->label &&
	       query_list[q]->n2->inGroup == g2->label) ||
	      (move_in == 2 &&
	       query_list[q]->n1->inGroup == g2->label &&
	       query_list[q]->n2->inGroup == oldg->label)) {
	    r--;
	    l -= IsThereLink(query_list[q]->n1, query_list[q]->n2);
	  }
	}
	dH += log(r + 1) + LogChoose(r, l);
	/* new configuration, new group */
	r = newg->size * g2->size;
	l = (*G2G)[newgnum][g2->label];
	for (q=0; q<nquery; q++) {
	  if ((move_in == 1 &&
	       query_list[q]->n1->inGroup == newg->label &&
	       query_list[q]->n2->inGroup == g2->label) ||
	      (move_in == 2 &&
	       query_list[q]->n1->inGroup == g2->label &&
	       query_list[q]->n2->inGroup == newg->label)) {
	    r--;
	    l -= IsThereLink(query_list[q]->n1, query_list[q]->n2);
	  }
	}
	dH += log(r + 1) + LogChoose(r, l);
      }
    }
    
    /* Metropolis acceptance */
    if ((dH <= 0.0) || (prng_get_next(gen) < exp(-dH))) {
      /* accept move: update energy */
      *H += dH;
    }
    else { 
      /* undo the move */
      MoveNode(node, newg, oldg);
      for (i=0; i<ngroup; i++) {
	(*G2G)[oldgnum][i] += n2gList[i];
	(*G2G)[newgnum][i] -= n2gList[i];
	(*G2Ginv)[i][oldgnum] += n2gList[i];
	(*G2Ginv)[i][newgnum] -= n2gList[i];
      }
    }
  } /* Moves completed: done! */

  return;
}


/*
  ---------------------------------------------------------------------
  Get the decorrelation step necessary to sample decorrelated
  partitions. From part1 and part2, the longest decorrelation step is
  chosen.
  ---------------------------------------------------------------------
*/
int
GetDecorrelationStep2State(double *H,
			   struct query **query_list, int nquery, 
			   struct node_gra **nlist1, struct node_gra **nlist2,
			   struct group **glist1, struct group **glist2,
			   struct group *part1, struct group *part2,
			   int nnod1, int nnod2,
			   int **G1G2, int **G2G1,
			   int *n2gList,
			   double **LogChooseList,
			   int LogChooseListSize,
			   struct prng *gen,
			   char verbose_sw)
{
  struct group *part1Ref, *part2Ref;
  int step, x1, x2;
  double y11, y12, y21, y22;
  double mutualInfo;
  int rep, nrep=10;
  double *decay1, meanDecay1, *decay2, meanDecay2;
  double **decay, meanDecay, sigmaDecay, result;
  int norm=0;

  x2 = (nnod1 + nnod2) / 5;
  if (x2 < 10)
    x2 = 10;
  x1 = x2 / 4;

  /* Get the nrep initial estimates */
  decay1 = allocate_d_vec(nrep);
  decay2 = allocate_d_vec(nrep);
  for (rep=0; rep<nrep; rep++) {
    switch (verbose_sw) {
    case 'q':
      break;
    default:
      fprintf(stderr, "#\n# Estimating decorrelation time (%d/%d)\n",
	      rep + 1, nrep);
      break;
    }
    part1Ref = CopyPartition(part1);
    part2Ref = CopyPartition(part2);
    for (step=0; step<=x2; step++) {
      MCStep2State(1, H, query_list, nquery, nlist1, nlist2,
		   glist1, glist2, part1, part2, nnod1, nnod2,
		   G1G2, G2G1, n2gList, LogChooseList, LogChooseListSize,
		   gen);
      if (step == x1)
	y11 = MutualInformation(part1Ref, part1);
	y12 = MutualInformation(part2Ref, part2);
    }
    y21 = MutualInformation(part1Ref, part1);
    y22 = MutualInformation(part2Ref, part2);
    if (nnod1 > 1)
      decay1[rep] = 2. * CalculateDecay(nnod1, x1, y11, x2, y21);
    else
      decay1[rep] = 1.e-6;
    if (nnod2 > 1)
      decay2[rep] = 2. * CalculateDecay(nnod2, x1, y12, x2, y22);
    else
      decay2[rep] = 1.e-6;
    switch (verbose_sw) {
    case 'q':
      break;
    default:
      fprintf(stderr, "# Decorrelation times (estimate %d) = %g %g\n",
	      rep + 1, decay1[rep], decay2[rep]);
      break;
    }
    if (decay1[rep] < 0. || decay2[rep] < 0.) {
      rep--;
      switch (verbose_sw) {
      case 'q':
	break;
      default:
	fprintf(stderr, "#\tignoring...\n");
	break;
      }
    }
    /* Free memory */
    RemovePartition(part1Ref);
    RemovePartition(part2Ref);
  }
  
  /* Get rid of bad estimates (Chauvenet criterion)  */
  meanDecay1 = mean(decay1, nrep);
  meanDecay2 = mean(decay2, nrep);
  if (meanDecay1 > meanDecay2) {
    meanDecay = meanDecay1;
    sigmaDecay = stddev(decay1, nrep);
    decay = &decay1;
    switch (verbose_sw) {
    case 'q':
      break;
    default:
      fprintf(stderr, "# Considering decorrelation from partition 1\n");
      break;
    }
  }
  else {
    meanDecay = meanDecay2;
    sigmaDecay = stddev(decay2, nrep);
    decay = &decay2;
    switch (verbose_sw) {
    case 'q':
      break;
    default:
      fprintf(stderr, "# Considering decorrelation from partition 2\n");
      break;
    }
  }
  result = meanDecay * nrep;
  for (rep=0; rep<nrep; rep++) {
    if (fabs((*decay)[rep] - meanDecay) / sigmaDecay > 2) {
      result -= (*decay)[rep];
      switch (verbose_sw) {
      case 'q':
	break;
      default:
	fprintf(stderr, "# Disregarding estimate %d\n", rep + 1);
	break;
      }
    }
    else {
      norm++;
    }
  }
  
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "# Decorrelation step: %d\n", (int)(result / norm + 0.5));
    break;
  }

  /* Clean up */
  free_d_vec(decay1);
  free_d_vec(decay2);

  return (int)(result / norm + 0.5);
}

/*
  ---------------------------------------------------------------------
  
  ---------------------------------------------------------------------
*/
void
ThermalizeMC2State(int decorStep,
		   double *H,
		   struct query **query_list, int nquery, 
		   struct node_gra **nlist1, struct node_gra **nlist2,
		   struct group **glist1, struct group **glist2,
		   struct group *part1, struct group *part2,
		   int nnod1, int nnod2,
		   int **G1G2, int **G2G1,
		   int *n2gList,
		   double **LogChooseList,
		   int LogChooseListSize,
		   struct prng *gen,
		   char verbose_sw)
{
  double HMean0=1.e10, HStd0=1.e-10, HMean1, HStd1, *Hvalues;
  int rep, nrep=20;
  int equilibrated=0;

  Hvalues = allocate_d_vec(nrep);

  do {
    
    /* MC steps */
    for (rep=0; rep<nrep; rep++) {
      MCStep2State(decorStep, H, query_list, nquery, nlist1, nlist2,
		   glist1, glist2, part1, part2, nnod1, nnod2,
		   G1G2, G2G1, n2gList, LogChooseList, LogChooseListSize,
		   gen);
      switch (verbose_sw) {
      case 'q':
	break;
      default:
	fprintf(stderr, "%lf\n", *H);
	break;
      }
      Hvalues[rep] = *H;
    }

    /* Check for equilibration */
    HMean1 = mean(Hvalues, nrep);
    HStd1 = stddev(Hvalues, nrep);
    if (HMean0 - HStd0 / sqrt(nrep) < HMean1 + HStd1 / sqrt(nrep)) {
      equilibrated++;
      switch (verbose_sw) {
      case 'q':
	break;
      default:
	fprintf(stderr, "#\tequilibrated (%d/5) H=%lf\n",
		equilibrated, HMean1);
	break;
      }
    }
    else {
      switch (verbose_sw) {
      case 'q':
	break;
      default:
	fprintf(stderr, "#\tnot equilibrated yet H0=%g+-%g H1=%g+-%g\n",
		HMean0, HStd0 / sqrt(nrep), HMean1, HStd1 / sqrt(nrep));
	break;
      }
      HMean0 = HMean1;
      HStd0 = HStd1;
      equilibrated = 0;
    }

  } while (equilibrated < 5);
  
  /* Clean up */
  free_d_vec(Hvalues);

  return;
}

/*
  -----------------------------------------------------------------------------
  Return the score of a single unobserved link i-j, that is
  p(A_ij=1|A^O), for a 2-state system
  -----------------------------------------------------------------------------
*/
double
LinkScore2State(struct binet *binet,
		struct query *the_query,
		int nIter,
		struct prng *gen,
		char verbose_sw,
		int decorStep)
{
  int nnod1=CountNodes(binet->net1), nnod2=CountNodes(binet->net2);
  struct node_gra *net1=NULL, *net2=NULL;
  struct group *part1=NULL, *part2=NULL;
  struct node_gra *p1=NULL, *p2=NULL, *node=NULL;
  struct node_gra **nlist1=NULL, **nlist2=NULL;
  struct group **glist1=NULL, **glist2=NULL;
  struct group *lastg=NULL;
  double H;
  int iter;
  double score=0.0, Z=0.0;
  int i, j;
  int **G1G2=NULL, **G2G1=NULL;
  int *n2gList=NULL;
  int LogChooseListSize = 500;
  double **LogChooseList=InitializeFastLogChoose(LogChooseListSize);
  struct node_lis *n1=NULL, *n2=NULL;
  double weight, contrib;
  int dice;
  int r, l;
  double mutualInfo;
  struct query *query_list[1];
  int nquery, query_linked;

  /*
    PRELIMINARIES
  */
  /* The query */
  query_list[0] = the_query;
  nquery = 1;
  query_linked = IsThereLink(the_query->n1, the_query->n2);

  /* Map nodes and groups to a list for faster access */
  nlist1 = (struct node_gra **) calloc(nnod1, sizeof(struct node_gra *));
  glist1 = (struct group **) calloc(nnod1, sizeof(struct group *));
  lastg = part1 = CreateHeaderGroup();
  p1 = net1 = binet->net1;
  while ((p1 = p1->next) != NULL) {
    nlist1[p1->num] = p1;
    lastg = glist1[p1->num] = CreateGroup(lastg, p1->num);
  }
  fprintf(stderr, "HERE %d %d\n", nnod1, nnod2);
  nlist2 = (struct node_gra **) calloc(nnod2, sizeof(struct node_gra *));
  fprintf(stderr, "HERE %d %d\n", nnod1, nnod2);
  glist2 = (struct group **) calloc(nnod2, sizeof(struct group *));
  fprintf(stderr, "HERE %d %d\n", nnod1, nnod2);
  lastg = part2 = CreateHeaderGroup();
  p2 = net2 = binet->net2;
  while ((p2 = p2->next) != NULL) {
    nlist2[p2->num] = p2;
    lastg = glist2[p2->num] = CreateGroup(lastg, p2->num);
  }
  fprintf(stderr, "HERE %d %d\n", nnod1, nnod2);

 /* Place nodes in random partitions */
  p1 = net1;
  ResetNetGroup(net1);
  while ((p1 = p1->next) != NULL) {
    dice = floor(prng_get_next(gen) * (double)nnod1);
    AddNodeToGroup(glist1[dice], p1);
  }
  p2 = net2;
  ResetNetGroup(net2);
  while ((p2 = p2->next) != NULL) {
    dice = floor(prng_get_next(gen) * (double)nnod2);
    AddNodeToGroup(glist2[dice], p2);
  }

  /* Get the initial group-to-group links matrix */
  G1G2 = allocate_i_mat(nnod1, nnod2);
  G2G1 = allocate_i_mat(nnod2, nnod1);
  n2gList = allocate_i_vec(max(nnod1, nnod2)); /* This is only used in
						  subroutines, but
						  allocated here for
						  convenience */
  for (i=0; i<nnod1; i++) {
    for (j=0; j<nnod2; j++) {
      G1G2[i][j] = G2G1[j][i] = NG2GLinks(glist1[i], glist2[j]);
/*       fprintf(stderr, "G1G2[%d][%d]=%d\n", i+1, j+1, G1G2[i][j]); */
    }
  }

  /*
    GET READY FOR THE SAMPLING
  */
  H = H2State(part1, part2, query_list, nquery);

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
    decorStep = GetDecorrelationStep2State(&H, query_list, nquery,
					   nlist1, nlist2, glist1, glist2,
					   part1, part2, nnod1, nnod2,
					   G1G2, G2G1, n2gList,
					   LogChooseList, LogChooseListSize,
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
		     query_list, nquery, 
		     nlist1, nlist2, glist1, glist2,
		     part1, part2, nnod1, nnod2,
		     G1G2, G2G1, n2gList,
		     LogChooseList, LogChooseListSize, gen, verbose_sw);
  
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
    MCStep2State(decorStep, &H, query_list, nquery, nlist1, nlist2,
		 glist1, glist2, part1, part2, nnod1, nnod2,
		 G1G2, G2G1, n2gList, LogChooseList, LogChooseListSize,
		 gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf\n", iter, H, H2State(part1, part2,
						       query_list, nquery));
/*       FPrintPartition(stderr, part1, 0); */
/*       FPrintPartition(stderr, part2, 0); */
      break;
    }

    /* Update partition function */
    weight = exp(-H);
    Z += weight;

    /* Update the score */
    l = G1G2[the_query->n1->inGroup][the_query->n2->inGroup] - query_linked;
    r = glist1[the_query->n1->inGroup]->size * \
      glist2[the_query->n2->inGroup]->size - 1;
    contrib = weight * (float)(l + 1) / (float)(r + 2);
    score += contrib;

  }  /* End of iter loop */

  /* Normalize the score */
  score /= Z;

  /* Done */
  RemovePartition(part1);
  RemovePartition(part2);
  free(glist1);
  free(glist2);
  free(nlist1);
  free(nlist2);
  free_i_mat(G1G2, nnod1);
  free_i_mat(G2G1, nnod2);
  free_i_vec(n2gList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  
  return score;
}
