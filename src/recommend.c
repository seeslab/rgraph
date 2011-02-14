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
#include <gsl/gsl_rng.h>

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
  -----------------------------------------------------------------------------
  Given a (bipartite) network of ratings, return the number of
  unobserved pairs.
  -----------------------------------------------------------------------------
*/
int
CountUnobserved(struct binet *ratings)
{
  return (CountNodes(ratings->net1) * CountNodes(ratings->net2)) -
    NLinksBipart(ratings);
}

/*
  -----------------------------------------------------------------------------
  Given a (bipartite) network of ratings, return a list of unobserved
  pairs.
  -----------------------------------------------------------------------------
*/
struct query **
BuildUnobservedSet(struct binet *ratings)
{
  struct node_gra *p, *m;
  int nunobserved=CountUnobserved(ratings), n=0;
  struct query **unobservedSet=NULL;

  unobservedSet = (struct query **) calloc(nunobserved, sizeof(struct query *));
  p = ratings->net1;
  while ((p = p->next) != NULL) {
    m = ratings->net2;
    while ((m = m->next) != NULL) {
      if (IsThereLink(p, m) == 0) {
	unobservedSet[n++] = CreateQuery(p, m);
      }
    }
  }
  return unobservedSet;
}

/*
  -----------------------------------------------------------------------------
  Given a (bipartite) network of ratings, remove all ratings whose
  value is r.
  -----------------------------------------------------------------------------
*/
void
RemoveRatings(struct binet *ratings, int r)
{
  struct node_gra *p=NULL;
  struct node_lis *m=NULL;

  p = ratings->net1;
  while ((p = p->next) != NULL) {
    m = p->neig;
    while (m->next!= NULL) {
      if ((m->next)->weight == r) {
	RemoveLink(p, (m->next)->ref, 1);
      }
      else {
	m = m->next;
      }
    }
  }
  return;
}

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  I/O functions
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
/*
  -----------------------------------------------------------------------------
  From a file with format:

  p1 m1 0
  p1 m2 5
  ...

  (representing person 1 rating movie 1 with a 0, movie 2 with a 5,
  and so on) return a bipartite network with link weights representing
  the ratings.
  -----------------------------------------------------------------------------
*/
struct binet *
ReadRecommenderObservations(FILE *inFile)
{
  return FBuildNetworkBipart(inFile, 1, 0);
}

/*
  -----------------------------------------------------------------------------
  Read a file with a query list with format:

  p1 m1
  p1 m2
  p2 m3
  ...

  representing the pairs (person, movie) that we want to predict. If a
  node in a query is missing from the bipartite ratings network, it is
  added with no links.
  
  -----------------------------------------------------------------------------
*/
struct query **
ReadQueries(FILE *inFile, int nQueries, struct binet *binet)
{
  void *dict1, *dict2;
  char node1[MAX_LABEL_LENGTH], node2[MAX_LABEL_LENGTH];
  struct query **querySet;
  struct query *q=NULL;
  int nq;
  struct node_gra *n1=NULL, *n2=NULL;

  /* Create the search dictionaries */
  dict1 = MakeLabelDict(binet->net1);
  dict2 = MakeLabelDict(binet->net2);

  /* Allocate space for the queries */
  querySet = (struct query **) calloc(nQueries, sizeof(struct query *));

  /* Go through the file and create the queries */
  for (nq=0; nq<nQueries; nq++) {
    fscanf(inFile, "%s %s\n", &node1[0], &node2[0]);
    n1 = GetNodeDict(node1, dict1);
    /* add to binet if absent */
    if (n1 == NULL) {
      n1 = CreateNodeGraph(binet->net1, node1);
      fprintf(stderr,
	      "WARNING! query node %s in net1 was never observed: adding!\n",
	      n1->label);
      FreeLabelDict(dict1);
      dict1 = MakeLabelDict(binet->net1);
    }
    n2 = GetNodeDict(node2, dict2);
    /* add to binet if absent */
    if (n2 == NULL) {
      n2 = CreateNodeGraph(binet->net2, node2);
      fprintf(stderr,
	      "WARNING! query node %s in net2 was never observed: adding!\n",
	      n2->label);
      FreeLabelDict(dict2);
      dict2 = MakeLabelDict(binet->net2);
    }
    querySet[nq] = CreateQuery(n1, n2);
  }

  /* Free memory */
  FreeLabelDict(dict1);
  FreeLabelDict(dict2);

  /* Done */
  return querySet;
}


/*
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
  Recommender functions
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------

/*
  -----------------------------------------------------------------------------
  Hamiltoninan H for 2-state recommender. ignore_list contains the list
  of unobserved links
  -----------------------------------------------------------------------------
*/
double
H2State(struct group *part1, struct group *part2,
	struct query **ignore_list, int nignore)
{
  struct group *g1=NULL, *g2=NULL;
  double r, l, H=0.0;
  int q;
  int n1=NGroups(part1), n2=NGroups(part2), nn1, nn2;
  int **rig=NULL, **lig=NULL;
  
  /* Discounted (unobserved or ignored) links for each group pair */
  rig = allocate_i_mat(n1, n2);
  lig = allocate_i_mat(n1, n2);
  for (nn1=0; nn1<n1; nn1++) {
    for (nn2=0; nn2<n2; nn2++) {
      rig[nn1][nn2] = lig[nn1][nn2] = 0;
    }
  }
  for (q=0; q<nignore; q++) {
    rig[ignore_list[q]->n1->inGroup][ignore_list[q]->n2->inGroup] += 1;
    lig[ignore_list[q]->n1->inGroup][ignore_list[q]->n2->inGroup] +=
      IsThereLink(ignore_list[q]->n1, ignore_list[q]->n2);
  }
  
  g1 = part1;
  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
      g2 = part2;
      while ((g2 = g2->next) != NULL) {
	if (g2->size > 0) {
	  r = g1->size * g2->size - rig[g1->label][g2->label];
	  l = NG2GLinks(g1, g2) - lig[g1->label][g2->label];
	  H += log(r + 1) + LogChoose(r, l);
	}
      }
    }
  }

  /* Free memory and return */
  free_i_mat(rig, n1);
  free_i_mat(lig, n1);
  return H;
}


/*
  -----------------------------------------------------------------------------
  Do a Monte Carlo step for the 2-state recommender system (faster
  than earlier versions by calculating rig and lig only twice).
  -----------------------------------------------------------------------------
*/
void
MCStep2State(int factor,
	     double *H,
	     struct query **ignore_list, int nignore, 
	     struct node_gra **nlist1, struct node_gra **nlist2,
	     struct group **glist1, struct group **glist2,
	     struct group *part1, struct group *part2,
	     int nnod1, int nnod2,
	     int **G1G2, int **G2G1,
	     int *n2gList,
	     double **LogChooseList,
	     int LogChooseListSize,
	     gsl_rng *gen)
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
  int **rig, **lig, nn1, nn2;

  /* Ratio of moves in each of the sets */
  set_ratio = (double)(nnod1*nnod1-1) / (double)(nnod1*nnod1+nnod2*nnod2-2);

  /* Discounted (unobserved or ignored) links for each group pair */
  rig = allocate_i_mat(nnod1, nnod2);
  lig = allocate_i_mat(nnod1, nnod2);

  /* The moves */
  for (move=0; move<(nnod1+nnod2)*factor; move++) {
    /* The move */
    if (gsl_rng_uniform(gen) < set_ratio) { /* move in first set */
      move_in = 1;
      G2G = &G1G2;
      G2Ginv = &G2G1;
      nnod = nnod1;
      ngroup = nnod2;
      dice = floor(gsl_rng_uniform(gen) * (double)nnod1);
      node = nlist1[dice];
      oldgnum = node->inGroup;
      do {
	newgnum = floor(gsl_rng_uniform(gen) * (double)nnod1);
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
      dice = floor(gsl_rng_uniform(gen) * (double)nnod2);
      node = nlist2[dice];
      oldgnum = node->inGroup;
      do {
	newgnum = floor(gsl_rng_uniform(gen) * (double)nnod2);
      } while (newgnum == oldgnum);
      oldg = glist2[oldgnum];
      newg = glist2[newgnum];
      g = g2 = part1;
    }

    /* THE CHANGE OF ENERGY */
    /* OLD CONFIGURATION CONTRIBUTION */

    /* discounted (unobserved or ignored) links for each group pair */
    for (nn1=0; nn1<nnod1; nn1++)
      for (nn2=0; nn2<nnod2; nn2++)
	rig[nn1][nn2] = lig[nn1][nn2] = 0;
    for (q=0; q<nignore; q++) {
      rig[ignore_list[q]->n1->inGroup][ignore_list[q]->n2->inGroup] += 1;
      lig[ignore_list[q]->n1->inGroup][ignore_list[q]->n2->inGroup] +=
	IsThereLink(ignore_list[q]->n1, ignore_list[q]->n2);
    }
    
    dH = 0.0;
    while ((g=g->next) != NULL) {
      if (g->size > 0) {  /* group is not empty */
      	n2gList[g->label] = NLinksToGroup(node, g);
	/* old configuration, old group */
	if (move_in == 1) {
	  r = oldg->size * g->size - rig[oldgnum][g->label];
	  l = (*G2G)[oldgnum][g->label]  - lig[oldgnum][g->label];
	}
	else {
	  r = oldg->size * g->size - rig[g->label][oldgnum];
	  l = (*G2G)[oldgnum][g->label]  - lig[g->label][oldgnum];
	}
	dH -= log(r + 1) + LogChoose(r, l);
	/* old configuration, new group */
	if (move_in == 1) {
	  r = newg->size * g->size - rig[newgnum][g->label];
	  l = (*G2G)[newgnum][g->label] - lig[newgnum][g->label];
	}
	else {
	  r = newg->size * g->size - rig[g->label][newgnum];
	  l = (*G2G)[newgnum][g->label] - lig[g->label][newgnum];
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

    /* NEW CONFIGURATION CONTRIBUTION */

    /* discounted (unobserved or ignored) links for each group pair */
    for (nn1=0; nn1<nnod1; nn1++)
      for (nn2=0; nn2<nnod2; nn2++)
	rig[nn1][nn2] = lig[nn1][nn2] = 0;
    for (q=0; q<nignore; q++) {
      rig[ignore_list[q]->n1->inGroup][ignore_list[q]->n2->inGroup] += 1;
      lig[ignore_list[q]->n1->inGroup][ignore_list[q]->n2->inGroup] +=
	IsThereLink(ignore_list[q]->n1, ignore_list[q]->n2);
    }

    while ((g2=g2->next) != NULL) {
      if (g2->size > 0) {  /* group is not empty */
	/* new configuration, old group */
	if (move_in == 1) {
	  r = oldg->size * g2->size - rig[oldgnum][g2->label];
	  l = (*G2G)[oldgnum][g2->label]  - lig[oldgnum][g2->label];
	}
	else {
	  r = oldg->size * g2->size - rig[g2->label][oldgnum];
	  l = (*G2G)[oldgnum][g2->label]  - lig[g2->label][oldgnum];
	}
	dH += log(r + 1) + LogChoose(r, l);
	/* new configuration, new group */
	if (move_in == 1) {
	  r = newg->size * g2->size - rig[newgnum][g2->label];
	  l = (*G2G)[newgnum][g2->label]  - lig[newgnum][g2->label];
	}
	else {
	  r = newg->size * g2->size - rig[g2->label][newgnum];
	  l = (*G2G)[newgnum][g2->label]  - lig[g2->label][newgnum];
	}
	dH += log(r + 1) + LogChoose(r, l);
      }
    }
    
    /* Metropolis acceptance */
    if ((dH <= 0.0) || (gsl_rng_uniform(gen) < exp(-dH))) {
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

  free_i_mat(rig, nnod1);
  free_i_mat(lig, nnod1);
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
			   struct query **ignore_list, int nignore, 
			   struct node_gra **nlist1, struct node_gra **nlist2,
			   struct group **glist1, struct group **glist2,
			   struct group *part1, struct group *part2,
			   int nnod1, int nnod2,
			   int **G1G2, int **G2G1,
			   int *n2gList,
			   double **LogChooseList,
			   int LogChooseListSize,
			   gsl_rng *gen,
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
      fprintf(stderr, "# %d / %d\n", step, x2);
      MCStep2State(1, H, ignore_list, nignore, nlist1, nlist2,
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
		   struct query **ignore_list, int nignore, 
		   struct node_gra **nlist1, struct node_gra **nlist2,
		   struct group **glist1, struct group **glist2,
		   struct group *part1, struct group *part2,
		   int nnod1, int nnod2,
		   int **G1G2, int **G2G1,
		   int *n2gList,
		   double **LogChooseList,
		   int LogChooseListSize,
		   gsl_rng *gen,
		   char verbose_sw)
{
  double HMean0=1.e10, HStd0=1.e-10, HMean1, HStd1, *Hvalues;
  int rep, nrep=20;
  int equilibrated=0;

  Hvalues = allocate_d_vec(nrep);

  do {
    
    /* MC steps */
    for (rep=0; rep<nrep; rep++) {
      MCStep2State(decorStep, H, ignore_list, nignore, nlist1, nlist2,
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
  p(A_ij=1|A^O), for a 2-state system. All links other than i-j are
  supposed to be observed, so links in the bipartite network indicate
  links with value 1, and non-links in the bipartite network indicate
  links with value 0.
  -----------------------------------------------------------------------------
*/
double
LinkScore2State(struct binet *binet,
		struct query *the_query,
		int nIter,
		gsl_rng *gen,
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
  nlist2 = (struct node_gra **) calloc(nnod2, sizeof(struct node_gra *));
  glist2 = (struct group **) calloc(nnod2, sizeof(struct group *));
  lastg = part2 = CreateHeaderGroup();
  p2 = net2 = binet->net2;
  while ((p2 = p2->next) != NULL) {
    nlist2[p2->num] = p2;
    lastg = glist2[p2->num] = CreateGroup(lastg, p2->num);
  }

 /* Place nodes in random partitions */
  p1 = net1;
  ResetNetGroup(net1);
  while ((p1 = p1->next) != NULL) {
    dice = floor(gsl_rng_uniform(gen) * (double)nnod1);
    AddNodeToGroup(glist1[dice], p1);
  }
  p2 = net2;
  ResetNetGroup(net2);
  while ((p2 = p2->next) != NULL) {
    dice = floor(gsl_rng_uniform(gen) * (double)nnod2);
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


/*
  -----------------------------------------------------------------------------
  Return 1 if q is in set, and 0 otherwise
  -----------------------------------------------------------------------------
*/
int
IsQueryInSet(struct query *q, struct query **set, int nset)
{
  int i;

  for (i=0; i<nset; i++)
    if ((q->n1 == set[i]->n1) && (q->n2 == set[i]->n2))
      return 1;
  return 0;
}


/*
  -----------------------------------------------------------------------------
  Return the score p(A_ij=1|A^O) of a collection {(i,j)} querySet of
  links, for a 2-state system. The ratings are a bipartite network
  with links (corresponding to observations) that have values 0 or 1.
  -----------------------------------------------------------------------------
*/
double *
MultiLinkScore2State(struct binet *ratings,
		     struct query **querySet, int nquery,
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
  double *score, Z=0.0;
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
  struct query **ignoreSet;
  int nignore, q, qq;
  int *ignore_linked;
  struct query **unobservedSet;
  int nunobserved;
  struct binet *binet=NULL; 
  void *dict1=NULL, *dict2=NULL;
  int **rig=NULL, **lig=NULL;
  
  /*
    PRELIMINARIES
  */
  /* Build the set of unobserved pairs */
  fprintf(stderr, ">> Building the set of unobserved pairs...\n");
  nunobserved = CountUnobserved(ratings);
  unobservedSet = BuildUnobservedSet(ratings);

  /* Create a network containing only 1-ratings (that is, without the
     0 ratings), and MAP THE QUERYSET AND THE UNOBSERVEDSET TO THE NEW
     NETWORK */
  fprintf(stderr, ">> Creating execution binet...\n");
  binet = CopyBipart(ratings);
  RemoveRatings(binet, 0);
  dict1 = MakeLabelDict(binet->net1);
  dict2 = MakeLabelDict(binet->net2);
  for (q=0; q<nquery; q++) {
    querySet[q]->n1 = GetNodeDict(querySet[q]->n1->label, dict1);
    querySet[q]->n2 = GetNodeDict(querySet[q]->n2->label, dict2);
  }
  for (q=0; q<nunobserved; q++) {
    unobservedSet[q]->n1 = GetNodeDict(unobservedSet[q]->n1->label, dict1);
    unobservedSet[q]->n2 = GetNodeDict(unobservedSet[q]->n2->label, dict2);
  }

  /* Initialize scores */
  score = allocate_d_vec(nquery);
  for (q=0; q<nquery; q++)
    score[q] = 0.0;

  /* The ignore set (we allocate the maximum possibly needed space) */
  fprintf(stderr, ">> Building the ignored set...\n");
  ignoreSet = (struct query **) calloc(nunobserved + nquery,
				       sizeof(struct query *));
  ignore_linked = allocate_i_vec(nunobserved + nquery);
  nignore = 0;
  for (q=0; q<nunobserved; q++) {
    ignoreSet[nignore] = unobservedSet[q];
    ignore_linked[nignore] = 0;
    nignore++;
  }
  for (q=0; q<nquery; q++) {
    if (IsQueryInSet(querySet[q], unobservedSet, nunobserved) == 0) {
      ignoreSet[nignore] = querySet[q];
      ignore_linked[nignore] = IsThereLink(querySet[q]->n1, querySet[q]->n2);
      nignore++;
    }
  }

  /* Map nodes and groups to a list for faster access */
  fprintf(stderr, ">> Mapping nodes and groups to lists...\n");
  nlist1 = (struct node_gra **) calloc(nnod1, sizeof(struct node_gra *));
  glist1 = (struct group **) calloc(nnod1, sizeof(struct group *));
  lastg = part1 = CreateHeaderGroup();
  p1 = net1 = binet->net1;
  while ((p1 = p1->next) != NULL) {
    nlist1[p1->num] = p1;
    lastg = glist1[p1->num] = CreateGroup(lastg, p1->num);
  }
  nlist2 = (struct node_gra **) calloc(nnod2, sizeof(struct node_gra *));
  glist2 = (struct group **) calloc(nnod2, sizeof(struct group *));
  lastg = part2 = CreateHeaderGroup();
  p2 = net2 = binet->net2;
  while ((p2 = p2->next) != NULL) {
    nlist2[p2->num] = p2;
    lastg = glist2[p2->num] = CreateGroup(lastg, p2->num);
  }

  /* Place nodes in random partitions */
  fprintf(stderr, ">> Placing nodes in random partitions...\n");
  p1 = net1;
  ResetNetGroup(net1);
  while ((p1 = p1->next) != NULL) {
    dice = floor(gsl_rng_uniform(gen) * (double)nnod1);
    AddNodeToGroup(glist1[dice], p1);
  }
  p2 = net2;
  ResetNetGroup(net2);
  while ((p2 = p2->next) != NULL) {
    dice = floor(gsl_rng_uniform(gen) * (double)nnod2);
    AddNodeToGroup(glist2[dice], p2);
  }

  /* Get the initial group-to-group links matrix */
  fprintf(stderr, ">> Getting the initial group-to-group links matrix...\n");
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

  /* Discounted (unobserved or ignored) links for each group pair */
  rig = allocate_i_mat(nnod1, nnod2);
  lig = allocate_i_mat(nnod1, nnod2);

  /*
    GET READY FOR THE SAMPLING
  */
  fprintf(stderr, ">> Getting the initial energy (nignore=%d)...\n", nignore);
  H = H2State(part1, part2, ignoreSet, nignore);

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
    decorStep = GetDecorrelationStep2State(&H, ignoreSet, nignore,
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
		     ignoreSet, nignore, 
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
    MCStep2State(decorStep, &H, ignoreSet, nignore, nlist1, nlist2,
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
						       ignoreSet, nignore));
      FPrintPartition(stderr, part1, 0);
      FPrintPartition(stderr, part2, 0);
      break;
    }

    /* Update partition function */
    weight = exp(-H);
    Z += weight;

    /* Update the scores */
    for (nn1=0; nn1<nnod1; nn1++) {
      for (nn2=0; nn2<nnod2; nn2++) {
	rig[nn1][nn2] = lig[nn1][nn2] = 0;
      }
    }
    for (q=0; q<nignore; q++) {
      rig[ignoreSet[q]->n1->inGroup][ignoreSet[q]->n2->inGroup] += 1;
      lig[ignoreSet[q]->n1->inGroup][ignoreSet[q]->n2->inGroup] +=
	IsThereLink(ignoreSet[q]->n1, ignoreSet[q]->n2);
    }
    for (q=0; q<nquery; q++) {
      l = G1G2[querySet[q]->n1->inGroup][querySet[q]->n2->inGroup] -
	lig[querySet[q]->n1->inGroup][querySet[q]->n2->inGroup];
      r = glist1[querySet[q]->n1->inGroup]->size *
	glist2[querySet[q]->n2->inGroup]->size -
	rig[querySet[q]->n1->inGroup][querySet[q]->n2->inGroup];
      contrib = weight * (float)(l + 1) / (float)(r + 2);
      score[q] += contrib;
    }

  }  /* End of iter loop */

  /* Normalize the scores */
  for (q=0; q<nquery; q++)
    score[q] /= Z;

  /*   Remap the queries to the original network */
  FreeLabelDict(dict1);  
  FreeLabelDict(dict2);  
  dict1 = MakeLabelDict(ratings->net1);
  dict2 = MakeLabelDict(ratings->net2);
  for (q=0; q<nquery; q++) {
    querySet[q]->n1 = GetNodeDict(querySet[q]->n1->label, dict1);
    querySet[q]->n2 = GetNodeDict(querySet[q]->n2->label, dict2);
  }
  FreeLabelDict(dict1);
  FreeLabelDict(dict2);

  /* Done */
  RemovePartition(part1);
  RemovePartition(part2);
  free(glist1);
  free(glist2);
  free(nlist1);
  free(nlist2);
  free_i_mat(G1G2, nnod1);
  free_i_mat(G2G1, nnod2);
  free_i_mat(rig, n1);
  free_i_mat(lig, n1);
  free_i_vec(n2gList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  free(ignoreSet);
  for (i=0; i<nunobserved; i++)
    FreeQuery(unobservedSet[i]);
  free(unobservedSet);
  free_i_vec(ignore_linked);
  RemoveBipart(binet);

  return score;
}
