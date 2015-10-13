/*
  missing.c
  $LastChangedDate: 2008-10-22 17:39:18 -0500 (Wed, 22 Oct 2008) $
  $Revision: 134 $
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf_gamma.h>

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "models.h"
#include "recommend.h"
#include "missing.h"


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
struct pair
{
  double x0;
  double y0;
  double x1;
  double y1;
};

int
ExponentialRootF(const gsl_vector *params,
		 void *points,
		 gsl_vector *f)
{
  const double a = gsl_vector_get(params, 0);
  const double b = gsl_vector_get(params, 1);

  double x0 = ((struct pair *)points)->x0;
  double y0 = ((struct pair *)points)->y0;
  double x1 = ((struct pair *)points)->x1;
  double y1 = ((struct pair *)points)->y1;
  
  const double r0 = y0 - a - (1. - a) * exp(-x0 / b);
  const double r1 = y1 - a - (1. - a) * exp(-x1 / b);
  
  gsl_vector_set(f, 0, r0);
  gsl_vector_set(f, 1, r1);
  
  return GSL_SUCCESS;
}

double
CalculateDecay(int nnod, double x1, double y1, double x2, double y2)
{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  int status;
  size_t i, iter = 0;
  const size_t n = 2;
  struct pair p = {x1, y1, x2, y2};
  gsl_multiroot_function f = {&ExponentialRootF, n, &p};
  double x_init[2] = {y2, sqrt(nnod)};
  gsl_vector *x = gsl_vector_alloc(n);
  double result;
  
  for (i=0; i<n; i++)
    gsl_vector_set(x, i, x_init[i]);
  
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);
  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      if (status)   /* check if solver is stuck */
	break;
      status =gsl_multiroot_test_residual(s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);
  
/*   fprintf(stderr, "# GSL status = %s\n", gsl_strerror(status)); */
  if (strcmp(gsl_strerror(status), "success") != 0)
    result = -1;
  else
    result = gsl_vector_get(s->x, 1);

  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return result;
}

/*
  ---------------------------------------------------------------------
  Partition H
  ---------------------------------------------------------------------
*/
double
PartitionH(struct group *part, double linC)
{
  struct group *g1=part, *g2;
  double r, l, H=0.0;
  int nnod=0, ng=0;

  H += linC * NNonEmptyGroups(part);

  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
      ng++;
      nnod += g1->size;
      r = g1->size * (g1->size - 1) / 2;
      l = g1->inlinks;
      H += log(r + 1) + LogChoose(r, l);
      g2 = g1;
      while ((g2 = g2->next) != NULL) {
	if (g2->size > 0) {
	  r = g1->size * g2->size;
	  l = NG2GLinks(g1, g2);
	  H += log(r + 1) + LogChoose(r, l);
	}
      }
    }
  }

  H -= gsl_sf_lnfact(nnod - ng);

  return H;
}

/*
  ---------------------------------------------------------------------
  Do a Monte Carlo step for the prediction of missing links. Each
  "step" involves nnod independent attempts to move a single node.
  ---------------------------------------------------------------------
*/
void
LinkScoreMCStep(int factor,
		double *H,
		double linC,
		struct node_gra **nlist,
		struct group **glist,
		struct group *part,
		int nnod,
		int **G2G,
		int *n2gList,
		double **LogChooseList,
		int LogChooseListSize,
		double *LogFactList, int LogFactListSize,
		gsl_rng *gen)
{
  double dH;
  struct group *g, *oldg, *newg;
  int dice, oldgnum, newgnum, n2g, oldg2g, newg2g, ng, noldg, nnewg;
  struct node_gra *node;
  int n2oldg, n2newg, newg2oldg, oldg2oldg, newg2newg;
  int r, l;
  int i, j, ngroup=nnod;
  int move;
  int effectng;
  /* Shortlist of groups. Groups not in this list are empty for
     sure. Note that, during a step, we may have to add groups to this
     list (groups that were initially empty but got filled) but we
     never remove groups from it (even if they loose all their nodes
     during a step). slgXsize is the size of slgX. */
  int *slg=NULL;
  int slgn, slgsize, isinlist;

  /* Preliminaries */
  slg = allocate_i_vec(nnod);  // Shortlist of groups
  slgsize = 0;
  g = part;
  while ((g=g->next) != NULL)
    if (g->size > 0)
      slg[slgsize++] = g->label;

  /* Steps */
  for (move=0; move<nnod*factor; move++) {

    /* Choose node and destination group */
    dice = floor(gsl_rng_uniform(gen) * (double)nnod);
    node = nlist[dice];
    oldgnum = node->inGroup;
    do {
      newgnum = floor(gsl_rng_uniform(gen) * (double)ngroup);
    } while (newgnum == oldgnum);
    oldg = glist[oldgnum];
    newg = glist[newgnum];
    
    /* Calculate the change of energy */
    dH = 0.0;
    noldg = oldg->size;
    nnewg = newg->size;
    if (noldg == 1)  /* number of groups would decrease by one */ 
      dH -= linC;
    if (nnewg == 0)  /* number of groups would increase by one */
      dH += linC;
    n2oldg = NLinksToGroup(node, oldg);
    n2newg = NLinksToGroup(node, newg);
    newg2oldg = NG2GLinks(newg, oldg);
    oldg2oldg = oldg->inlinks;
    newg2newg = newg->inlinks;
    effectng = 0;
    for (slgn=0; slgn<slgsize; slgn++) {
      g = glist[slg[slgn]];
      if (g->size > 0) {  /* group is not empty */
	effectng++;
	n2gList[g->label] = NLinksToGroup(node, g);
	if (g->label == oldg->label) {
	  /* old conf, oldg-oldg */
	  r = noldg * (noldg - 1) / 2;
	  l = oldg2oldg;
	  dH -= log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /* new conf, oldg-olg */
	  r = (noldg - 1) * (noldg - 2) / 2;
	  l = oldg2oldg - n2oldg;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /* old conf, newg-oldg */
	  r = nnewg * noldg;
	  l = newg2oldg;
	  dH -= log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /* new conf, newg-oldg */
	  r = (nnewg + 1) * (noldg - 1);
	  l = newg2oldg + n2oldg - n2newg;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	}
	else if (g->label == newg->label) {
	  /* old conf, newg-newg */
	  r = nnewg * (nnewg - 1) / 2;
	  l = newg2newg;
	  dH -= log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /* new conf, newg-olg */
	  r = (nnewg + 1) * nnewg / 2;
	  l = newg2newg + n2newg;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	}
	else {
	  n2g = n2gList[g->label];
	  oldg2g = G2G[oldg->label][g->label];
	  newg2g = G2G[newg->label][g->label];
	  ng = g->size;
	  /* old conf, oldg-g */
	  r = noldg * ng;
	  l = oldg2g;
	  dH -= log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /* new conf, oldg-g */
	  r = (noldg - 1) * ng;
	  l = oldg2g - n2g;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /* old conf, newg-g */
	  r = nnewg * ng;
	  l = newg2g;
	  dH -= log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /* new conf, newg-g */
	  r = (nnewg + 1) * ng;
	  l = newg2g + n2g;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	}
      }
      else { /* group is empty */
	n2gList[g->label] = 0.0;
      }
    }

    /* Labeled-groups sampling correction */
    if (noldg == 1 && nnewg != 0) { // effectng would decrease by one
      dH -= -FastLogFact(ngroup - effectng, LogFactList, LogFactListSize);
      dH += -FastLogFact(ngroup - (effectng-1), LogFactList, LogFactListSize);
    }
    else if (noldg != 1 && nnewg == 0) {  // effectng would increase by one
      dH -= -FastLogFact(ngroup - effectng, LogFactList, LogFactListSize);
      dH += -FastLogFact(ngroup - (effectng+1), LogFactList, LogFactListSize);
    }

    /* Metropolis rule */
    if ((dH <= 0.0) || (gsl_rng_uniform(gen) < exp(-dH))) {
      
      /* accept move */
      MoveNode(node, oldg, newg);
      *H += dH;
    
      /* update G2G */
      for (i=0; i<nnod; i++) {
	G2G[i][oldgnum] -= n2gList[i];
	G2G[oldgnum][i] = G2G[i][oldgnum];
	G2G[i][newgnum] += n2gList[i];
	G2G[newgnum][i] = G2G[i][newgnum];
      }

      /* Add newg to shortlist (if it's not there already!) */
      isinlist = 0;
      for (slgn=0; slgn<slgsize; slgn++)
	if (newgnum == slg[slgn])
	  isinlist = 1;
      if (isinlist == 0)
	slg[slgsize++] = newgnum;
    }
  } /* nnod moves completed: done! */

  free_i_vec(slg);
}


/*
  ---------------------------------------------------------------------
  
  ---------------------------------------------------------------------
*/
int
GetDecorrelationStep(double *H,
		     double linC,
		     struct node_gra **nlist,
		     struct group **glist,
		     struct group *part,
		     int nnod,
		     int **G2G,
		     int *n2gList,
		     double **LogChooseList,
		     int LogChooseListSize,
		     double *LogFactList, int LogFactListSize,
		     gsl_rng *gen,
		     char verbose_sw)
{
  struct group *partRef;
  int step, x1, x2;
  double y1, y2;
  double mutualInfo;
  int rep, nrep=10;
  double *decay, meanDecay, sigmaDecay, result;
  int norm=0;

  x2 = nnod / 5;
  x1 = x2 / 4;

  /* Get the nrep initial estimates */
  decay = allocate_d_vec(nrep);
  for (rep=0; rep<nrep; rep++) {
    switch (verbose_sw) {
    case 'q':
      break;
    default:
      fprintf(stderr, "#\n# Estimating decorrelation time (%d/%d)\n",
	      rep + 1, nrep);
      break;
    }
    partRef = CopyPartition(part);
    for (step=0; step<=x2; step++) {
      LinkScoreMCStep(1, H, linC, nlist, glist, part,
		      nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
		      LogFactList, LogFactListSize, 
		      gen);
      if (step == x1)
	y1 = MutualInformation(partRef, part, 0);
    }
    y2 = MutualInformation(partRef, part, 0);
    RemovePartition(partRef);
    decay[rep] = 2. * CalculateDecay(nnod, x1, y1, x2, y2);
    switch (verbose_sw) {
    case 'q':
      break;
    default:
      fprintf(stderr, "# Decorrelation time (estimate %d) = %g\n",
	      rep + 1, decay[rep]);
      break;
    }
    if (decay[rep] < 0) {
      rep--;
      switch (verbose_sw) {
      case 'q':
	break;
      default:
	fprintf(stderr, "#\tignoring...\n");
	break;
      }
    }
  }
  
  /* Get rid of bad estimates (Chauvenet criterion)  */
  meanDecay = mean(decay, nrep);
  sigmaDecay = stddev(decay, nrep);
  result = meanDecay * nrep;
  for (rep=0; rep<nrep; rep++) {
    if (fabs(decay[rep] - meanDecay) / sigmaDecay > 2) {
      result -= decay[rep];
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
  
  /* Clean up */
  free_d_vec(decay);

  return result / norm;
}

/*
  ---------------------------------------------------------------------
  
  ---------------------------------------------------------------------
*/
void
ThermalizeLinkScoreMC(int decorStep,
		      double *H,
		      double linC,
		      struct node_gra **nlist,
		      struct group **glist,
		      struct group *part,
		      int nnod,
		      int **G2G,
		      int *n2gList,
		      double **LogChooseList,
		      int LogChooseListSize,
		      double *LogFactList, int LogFactListSize,
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
      LinkScoreMCStep(decorStep, H, linC, nlist, glist, part,
		      nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
		      LogFactList, LogFactListSize,
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
    if ((HMean0 - HStd0 / sqrt(nrep)) - (HMean1 + HStd1 / sqrt(nrep))
	< EPSILON) {
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
  ---------------------------------------------------------------------
  Predict the links missing in a network. The algorithm returns a
  matrix of scores for all links.
  ---------------------------------------------------------------------
*/
double **
LinkScore(struct node_gra *net,
	  double linC,
	  int nIter,
	  gsl_rng *gen,
	  char verbose_sw)
{
  int nnod=CountNodes(net);
  struct group *part=NULL;
  struct node_gra *p=NULL, *node=NULL;
  struct node_gra **nlist=NULL;
  struct group **glist=NULL;
  struct group *lastg=NULL;
  double H;
  int iter, decorStep;
  double **predA=NULL;
  int i, j;
  int **G2G=NULL;
  int *n2gList=NULL;
  int LogChooseListSize = 500;
  double **LogChooseList=InitializeFastLogChoose(LogChooseListSize);
  struct node_lis *p1=NULL, *p2=NULL;
  double contrib;
  int dice;
  int norm = 0;
  int r, l;
  double mutualInfo;
  int LogFactListSize = 10000;
  double *LogFactList=InitializeFastLogFact(LogFactListSize);

  /*
    PRELIMINARIES
  */
  /* Initialize the predicted adjacency matrix */
  predA = allocate_d_mat(nnod, nnod);
  for (i=0; i<nnod; i++)
    for (j=0; j<nnod; j++)
      predA[i][j] = 0.0;

  /* Map nodes and groups to a list for faster access */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  glist = (struct group **) calloc(nnod, sizeof(struct group *));
  lastg = part = CreateHeaderGroup();
  p = net;
  while ((p = p->next) != NULL) {
    nlist[p->num] = p;
    lastg = glist[p->num] = CreateGroup(lastg, p->num);
  }

  /* Place nodes in random partitions */
  p = net;
  ResetNetGroup(net);
  while ((p = p->next) != NULL) {
    dice = floor(gsl_rng_uniform(gen) * (double)nnod);
    AddNodeToGroup(glist[dice], p);
  }

  /* Get the initial group-to-group links matrix */
  G2G = allocate_i_mat(nnod, nnod);
  n2gList = allocate_i_vec(nnod);
  for (i=0; i<nnod; i++) {
    G2G[i][i] = glist[i]->inlinks;
    for (j=i+1; j<nnod; j++) {
      G2G[i][j] = G2G[j][i] = NG2GLinks(glist[i], glist[j]);
    }
  }

  /*
    GET READY FOR THE SAMPLING
  */
  /* Get the decorrelation time */
  H = PartitionH(part, linC);
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "# CALCULATING DECORRELATION TIME\n");
    fprintf(stderr, "# ------------------------------\n");
    break;
  }
  decorStep = GetDecorrelationStep(&H, linC, nlist, glist, part,
				   nnod, G2G, n2gList,
				   LogChooseList, LogChooseListSize,
				   LogFactList, LogFactListSize,
				   gen, verbose_sw);
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "\n# Decorrelation step = %d\n\n", decorStep);
    break;
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
  ThermalizeLinkScoreMC(decorStep, &H, linC, nlist, glist, part,
			nnod, G2G, n2gList,
			LogChooseList, LogChooseListSize,
			LogFactList, LogFactListSize,
			gen, verbose_sw);
  
  /*
    SAMPLIN' ALONG
  */
  /* Unless we are in debug mode, reset the origin of energies */
  switch (verbose_sw) {
  case 'd':
    break;
  default:
    /* H = 0; */
    break;
  }

  /* Do the MC Steps */
  for (iter=0; iter<nIter; iter++) {
    LinkScoreMCStep(decorStep, &H, linC, nlist, glist, part,
		    nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
		    LogFactList, LogFactListSize,
		    gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf\n", iter, H, PartitionH(part, linC));
      FPrintPartition(stderr, part, 0);
      break;
    }

    /* Update partition function */
    norm += 1;

    /* Update the predicted adjacency matrix by going through all
       group pairs */
    for (i=0; i<nnod; i++) {
      if (glist[i]->size > 0) {
	
	/* update the within-group pairs */
	r = glist[i]->size * (glist[i]->size - 1) / 2;
	l = glist[i]->inlinks;
	contrib = (float)(l + 1) / (float)(r + 2);
	p1 = glist[i]->nodeList;
	while ((p1 = p1->next) != NULL) {
	  p2 = p1;
	  while ((p2 = p2->next) != NULL) {
	    predA[p1->node][p2->node] += contrib;
	    predA[p2->node][p1->node] += contrib;
	  }
	}
      
	/* update the between-group pairs */
	for (j=i+1; j<nnod; j++) {
	  if (glist[j]->size > 0) {
	    l = G2G[i][j];
	    r = glist[i]->size * glist[j]->size;
	    contrib = (float)(l + 1) / (float)(r + 2);
	    p1 = glist[i]->nodeList;
	    while ((p1 = p1->next) != NULL) {
	      p2 = glist[j]->nodeList;
	      while ((p2 = p2->next) != NULL) {
		predA[p1->node][p2->node] += contrib;
		predA[p2->node][p1->node] += contrib;
	      }
	    }
	  }
	}
      }
    } /* Done updating adjacency matrix */

  }  /* End of iter loop */

  /* Normalize the predicted adjacency matrix */
  for (i=0; i<nnod; i++) {
    for (j=0; j<nnod; j++) {
      predA[i][j] /= (double)norm;
    }
  }

  /* Done */
  RemovePartition(part);
  free(glist);
  free(nlist);
  free_i_mat(G2G, nnod);
  free_i_vec(n2gList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  FreeFastLogFact(LogFactList);

  return predA;
}

/*
  ---------------------------------------------------------------------
  Return error of the stochastic block model for a given network
  ---------------------------------------------------------------------
*/
double
SBMError(struct node_gra *net, gsl_rng *gen)
{
  struct node_gra **nlist;
  int i, j, nnod=CountNodes(net);
  struct node_gra *p;
  double score=0.0;
  double **pairScore;

  /* Get the link score */
  pairScore = LinkScore(net, 0.0, 10000, gen, 'q');

  /* Map nodes to a list for faster access */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  p = net;
  while ((p = p->next) != NULL) {
    nlist[p->num] = p;
  }

  /* Calculate the reliability score */
  for (i=0; i<nnod; i++) {
    for (j=i+1; j<nnod; j++) {
/*       fprintf(stderr, "%s %s %lf\n", */
/* 	      nlist[i]->label, nlist[j]->label, pairScore[i][j]); */
      if (IsThereLink(nlist[i], nlist[j]) == 1) {
	score += (1.0 - pairScore[i][j]) * (1.0 - pairScore[i][j]);
      }
      else {
	score += pairScore[i][j] * pairScore[i][j];
      }
    }
  }

  /* Normalize */
  score = sqrt(score / (double)(nnod * (nnod - 1) / 2));

  /* Done */
  free(nlist);
  free_d_mat(pairScore, nnod);
  return score;
}

/*
  ---------------------------------------------------------------------
  Return normalized error of the stochastic block model for a given
  network
  ---------------------------------------------------------------------
*/
double
SBMStructureScore(struct node_gra *net, int nrep, gsl_rng *gen)
{
  int rep;
  double scoreOri, *scoreRan, score;
  struct node_gra *randNet;

  /* Original network */
  scoreOri = SBMError(net, gen);

  /* Randomizations */
  scoreRan = allocate_d_vec(nrep);
  randNet = CopyNetwork(net);
  for (rep=0; rep<nrep; rep++) {
    RandomizeSymmetricNetwork(randNet, 100, gen);
    scoreRan[rep] = SBMError(randNet, gen);
  }

  /* The score */
  score = (mean(scoreRan, nrep) - scoreOri) / stddev(scoreRan, nrep);

  /* Done */
  free_d_vec(scoreRan);
  RemoveGraph(randNet);
  return score;
}

/*
  ---------------------------------------------------------------------
  Create a network from the SBM link scores: if the score is q_ij>0.5,
  then the network has a link between i an j, and vice versa.
  ---------------------------------------------------------------------
*/
struct node_gra *
NetFromSBMScores(struct node_gra *net, gsl_rng *gen)
{
  struct node_gra **nlist;
  int n1, n2, nnod=CountNodes(net);
  struct node_gra *p_new, *p, *net_new;
  double score=0.0;
  double **pairScore;

  /* Get the link score */
  pairScore = LinkScore(net, 0.0, 10000, gen, 'q');

  /* Create an empty network */
  net_new = EmptyGraph(nnod);

  /* Map nodes in the new network to a list for faster access, and
     rename them as in the original network */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  p = net;
  p_new = net_new;
  while ((p_new = p_new->next) != NULL) {
    p = p->next;
    strcpy(p_new->label, p->label);
    nlist[p_new->num] = p_new;
  }

  /* Add the links */
  for (n1=0; n1<nnod; n1++) {
    for (n2=n1+1; n2<nnod; n2++) {
/*       if (gsl_rng_uniform(gen) < pairScore[n1][n2]) { */
      if (pairScore[n1][n2] > 0.5) {
	AddAdjacency(nlist[n1], nlist[n2], 0, 0, 0, 0);
	AddAdjacency(nlist[n2], nlist[n1], 0, 0, 0, 0);
      }
    }
  }

  /* Done */
  free(nlist);
  free_d_mat(pairScore, nnod);
  return net_new;
}

/*
  ---------------------------------------------------------------------
  Given a "target" network A and an observed network A^O, returns the
  network reliability R_A of A divided by the network reliability of
  the observation A^O. The normalization is done to avoid as much as
  possible huge numbers.
  ---------------------------------------------------------------------
*/
double
NetworkScore(struct node_gra *netTar,
	     struct node_gra *netObs,
	     double linC,
	     int nIter,
	     gsl_rng *gen,
	     char verbose_sw)
{
  int nnod=CountNodes(netObs);
  struct group *part=NULL, *partCopy=NULL;
  struct node_gra *p=NULL, *node=NULL;
  struct node_gra **nlist=NULL;
  struct group **glist=NULL;
  struct group *g1Obs=NULL, *g2Obs=NULL, *g1Tar=NULL, *g2Tar=NULL;
  struct group *lastg=NULL;
  double H;
  int iter, decorStep;
  int norm=0;
  int **G2G=NULL;
  int *n2gList=NULL;
  int LogChooseListSize = 500;
  double **LogChooseList=InitializeFastLogChoose(LogChooseListSize);
  struct node_lis *p1=NULL, *p2=NULL;
  double scoreTar=0.0, scoreObs=0.0;
  double contribObs, contribTar, contribBase;
  int contribBase_sw=0;
  int i, j, dice;
  int r, lObs, lTar;
  double mutualInfo;
  int LogFactListSize = 10000;
  double *LogFactList=InitializeFastLogFact(LogFactListSize);

  /*
    PRELIMINARIES
  */
  /* Map nodes and groups to a list for faster access */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  glist = (struct group **) calloc(nnod, sizeof(struct group *));
  lastg = part = CreateHeaderGroup();
  p = netObs;
  while ((p = p->next) != NULL) {
    nlist[p->num] = p;
    lastg = glist[p->num] = CreateGroup(lastg, p->num);
  }

  /* Place nodes in random partitions */
  p = netObs;
  ResetNetGroup(netObs);
  while ((p = p->next) != NULL) {
    dice = floor(gsl_rng_uniform(gen) * (double)nnod);
    AddNodeToGroup(glist[dice], p);
  }

  /* Get the initial group-to-group links matrix */
  G2G = allocate_i_mat(nnod, nnod);
  n2gList = allocate_i_vec(nnod);
  for (i=0; i<nnod; i++) {
    G2G[i][i] = glist[i]->inlinks;
    for (j=i+1; j<nnod; j++) {
      G2G[i][j] = G2G[j][i] = NG2GLinks(glist[i], glist[j]);
    }
  }

  /*
    GET READY FOR THE SAMPLING
  */
  /* Get the decorrelation time */
  H = PartitionH(part, linC);
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "# CALCULATING DECORRELATION TIME\n");
    fprintf(stderr, "# ------------------------------\n");
    break;
  }
  decorStep = GetDecorrelationStep(&H, linC, nlist, glist, part,
				   nnod, G2G, n2gList,
				   LogChooseList, LogChooseListSize,
				   LogFactList, LogFactListSize, 
				   gen, verbose_sw);

  /* Thermalization */
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "#\n#\n# THERMALIZING\n");
    fprintf(stderr, "# ------------\n");
  }
  ThermalizeLinkScoreMC(decorStep, &H, linC, nlist, glist, part,
			nnod, G2G, n2gList,
			LogChooseList, LogChooseListSize,
			LogFactList, LogFactListSize,
			gen, verbose_sw);
  
  /*
    SAMPLIN' ALONG
  */
  /* Unless we are in debug mode, reset the origin of energies to
     avoid huge exponentials */
  switch (verbose_sw) {
  case 'd':
    break;
  default:
    H = 0;
    break;
  }
  for (iter=0; iter<nIter; iter++) {
    LinkScoreMCStep(decorStep, &H, linC, nlist, glist, part,
		    nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
		    LogFactList, LogFactListSize,
		    gen);
    switch (verbose_sw) {
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf\n", iter, H, PartitionH(part, linC));
      break;
    }

    /* Update the partition function */
    norm += 1;

    /* Copy the partition and map it to the target network */
    partCopy = CopyPartition(part);
    MapPartToNet(partCopy, netTar);
    MapPartToNet(part, netObs);

    /* Update the score */
    g1Obs = part;
    g1Tar = partCopy;
    contribTar = contribObs = 0.0;
    while ((g1Obs = g1Obs->next) != NULL) {
      g1Tar = g1Tar->next;
      if (g1Obs->size > 0) {
	/* update the within-group pairs */
	r = g1Obs->size * (g1Obs->size - 1) / 2;
	lObs = g1Obs->inlinks;
	lTar = g1Tar->inlinks;
	contribTar += (+log(r + 1) + LogChoose(r, lObs)
		       -log(2 * r + 1) - LogChoose(2 * r, lObs + lTar));
	contribObs += (+log(r + 1) + LogChoose(r, lObs)
		       -log(2 * r + 1) - LogChoose(2 * r, 2 * lObs));
     
	/* update the between-group pairs */
	g2Obs = g1Obs;
	g2Tar = g1Tar;
	while ((g2Obs = g2Obs->next) != NULL) {
	  g2Tar = g2Tar->next;
	  if (g2Obs->size > 0) {
	    r = g1Obs->size * g2Obs->size;
	    lObs = G2G[g1Obs->label][g2Obs->label];
	    lTar = NG2GLinks(g1Tar, g2Tar);
	    contribTar += (+log(r + 1) + LogChoose(r, lObs)
			   -log(2 * r + 1) - LogChoose(2 * r, lObs + lTar));
	    contribObs += (+log(r + 1) + LogChoose(r, lObs)
			   -log(2 * r + 1) - LogChoose(2 * r, 2 * lObs));
	  }
	}
      }
    } /* Done calculating the contribution of this partition */
    if (contribBase_sw == 0) {
      contribBase = contribObs;
      contribBase_sw = 1;
    }      
    scoreTar += exp(contribTar - contribBase);
    scoreObs += exp(contribObs - contribBase);

    RemovePartition(partCopy);
  }  /* End of iter loop */

  /* Normalize the scores */
  scoreTar /= (double)norm;
  scoreObs /= (double)norm;

  /* Done */
  RemovePartition(part);
  free(glist);
  free(nlist);
  free_i_mat(G2G, nnod);
  free_i_vec(n2gList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  
  return scoreTar / scoreObs;
}


/*
  ---------------------------------------------------------------------
  Reconstruct a network using a up-hill search in the NetworkScore
  space
  ---------------------------------------------------------------------
*/
struct node_gra *
NetReconstruct(struct node_gra *netObs,
	       gsl_rng *gen)
{
  struct node_gra *net, **nlist, *p, *n1Add, *n2Add, *n1Remove, *n2Remove;
  int i, n1, n2;
  double **linkScores;
  double scoreNew, scoreOld, maxLinkScore, minLinkScore;
  int nnod=CountNodes(netObs);
  int someChanged, someRejected;
  int thresReject = 5;

  /* Map nodes to list for fast access */
  net = CopyNetwork(netObs);
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  p = net;
  while ((p = p->next) != NULL)
    nlist[p->num] = p;

  /* Get the initial network score */
  scoreOld = NetworkScore(net, netObs, 0.0, 100, gen, 'q');

  /* ITERATIVE PROCEDURE */
  do {   /* Keep going until no changes occur */

    /* Get the (quick'n'dirty) link scores */
    linkScores = LinkScore(net, 0.0, 1000, gen, 'q');

    someChanged = 0;
    someRejected = 0;
    do {   /* Keep going until thresReject consecutive changes are
	      rejected */

      /* Get the largest missing link score and the smallest link
	 score */
      maxLinkScore = 0.0;
      minLinkScore = 1.0;
      for (n1=0; n1<nnod; n1++) {
	for (n2=n1+1; n2<nnod; n2++) {
	  if ((linkScores[n1][n2] > maxLinkScore) &&
	      (IsThereLink(nlist[n1], nlist[n2]) == 0)) {
	    n1Add = nlist[n1];
	    n2Add = nlist[n2];
	    maxLinkScore = linkScores[n1][n2];
	  }
	  if ((linkScores[n1][n2] < minLinkScore) &&
	      (IsThereLink(nlist[n1], nlist[n2]) == 1)) {
	    n1Remove = nlist[n1];
	    n2Remove = nlist[n2];
	    minLinkScore = linkScores[n1][n2];
	  }
	}
      }

      /* Add and remove links */
      AddAdjacency(n1Add, n2Add, 0, 0, 1.0, 2);
      AddAdjacency(n2Add, n1Add, 0, 0, 1.0, 2);
      RemoveLink(n1Remove, n2Remove, 1);

      /* Calculate new score and do the comparison */
      scoreNew = NetworkScore(net, netObs, 0.0, 100, gen, 'q');
      if (scoreNew > scoreOld) {
/* 	fprintf(stderr, "Accepting %g %g\n", scoreOld, scoreNew); */
	scoreOld = scoreNew;
	someRejected = 0;
	someChanged = 1;
      }
      else {
/* 	fprintf(stderr, "Rejecting %g %g\n", scoreOld, scoreNew); */
	someRejected += 1;
	/* undo the change */
	AddAdjacency(n1Remove, n2Remove, 0, 0, 1.0, 3);
	AddAdjacency(n2Remove, n1Remove, 0, 0, 1.0, 3);
	RemoveLink(n1Add, n2Add, 1);
      }
    } while(someRejected < thresReject);
  } while (someChanged == 1);

  return net;
}


/*
  ---------------------------------------------------------------------
  Return all partitions sampled according to the reliability formalism
  ---------------------------------------------------------------------
*/
struct group **
PartitionSampling(struct node_gra *net,
		  double linC,
		  int nIter,
		  gsl_rng *gen,
		  char verbose_sw,
		  int burnin,
		  double thinning)
{
  int nnod=CountNodes(net);
  struct group *part=NULL;
  struct node_gra *p=NULL, *node=NULL;
  struct node_gra **nlist=NULL;
  struct group **glist=NULL;
  struct group **partList=NULL;
  struct group *lastg=NULL;
  double H;
  int iter, decorStep;
  int **G2G=NULL;
  int *n2gList=NULL;
  int LogChooseListSize = 500;
  double **LogChooseList=InitializeFastLogChoose(LogChooseListSize);
  int LogFactListSize = 10000;
  double *LogFactList=InitializeFastLogFact(LogFactListSize);
  int dice;
  int i, j;

  /*
    PRELIMINARIES
  */
  /* Allocate memory for the list of partitions */
  partList = (struct group **) calloc(nIter, sizeof(struct group *));

  /* Map nodes and groups to a list for faster access */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  glist = (struct group **) calloc(nnod, sizeof(struct group *));
  lastg = part = CreateHeaderGroup();
  p = net;
  while ((p = p->next) != NULL) {
    nlist[p->num] = p;
    lastg = glist[p->num] = CreateGroup(lastg, p->num);
  }

  /* Place nodes in random partitions */
  p = net;
  ResetNetGroup(net);
  while ((p = p->next) != NULL) {
    dice = floor(gsl_rng_uniform(gen) * (double)nnod);
    AddNodeToGroup(glist[dice], p);
  }

  /* Get the initial group-to-group links matrix */
  G2G = allocate_i_mat(nnod, nnod);
  n2gList = allocate_i_vec(nnod);
  for (i=0; i<nnod; i++) {
    G2G[i][i] = glist[i]->inlinks;
    for (j=i+1; j<nnod; j++) {
      G2G[i][j] = G2G[j][i] = NG2GLinks(glist[i], glist[j]);
    }
  }

  /*
    GET READY FOR THE SAMPLING
  */
  /* Get the decorrelation time */
  H = PartitionH(part, linC);
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "# CALCULATING DECORRELATION TIME\n");
    fprintf(stderr, "# ------------------------------\n");
    break;
  }
  decorStep = GetDecorrelationStep(&H, linC, nlist, glist, part,
				   nnod, G2G, n2gList,
				   LogChooseList, LogChooseListSize,
				   LogFactList, LogFactListSize,
				   gen, verbose_sw);

  /* Thermalization and additional burn-in*/
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "#\n#\n# THERMALIZING\n");
    fprintf(stderr, "# ------------\n");
    break;
  }
  ThermalizeLinkScoreMC(decorStep, &H, linC, nlist, glist, part,
			nnod, G2G, n2gList,
			LogChooseList, LogChooseListSize,
			LogFactList, LogFactListSize,
			gen, verbose_sw);
  for (iter=0; iter<burnin; iter++) {
    LinkScoreMCStep(decorStep, &H, linC, nlist, glist, part,
		    nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
		    LogFactList, LogFactListSize,
		    gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "burn %d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "burn %d %lf %lf\n", iter, H, PartitionH(part, linC));
      break;
    }
  }
  
  /*
    SAMPLIN' ALONG
  */
  /* Unless we are in debug mode, reset the origin of energies */
  switch (verbose_sw) {
  case 'd':
    break;
  default:
    H = 0;
    break;
  }

  /* Do the MC Steps */
  for (iter=0; iter<nIter; iter++) {
    LinkScoreMCStep((int)(decorStep*thinning), &H, linC, nlist, glist, part,
		    nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
		    LogFactList, LogFactListSize,
		    gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf\n", iter, H, PartitionH(part, linC));
      break;
    }

    /* Save the partition */
    partList[iter] = CopyPartition(part);
    MapPartToNet(part, net);

  }  /* End of iter loop */

  /* Done */
  RemovePartition(part);
  free(glist);
  free(nlist);
  free_i_mat(G2G, nnod);
  free_i_vec(n2gList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  
  return partList;
}



/*
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
  Link reliability K-state
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
*/

/*
  -----------------------------------------------------------------------------
  Hamiltoninan H for K-state reliability estimation. The part
  partition must be mapped onto a network whose links represent
  observed ratings. The weight of the link represents the
  rating. CAUTION!! ALL LINKS IN THE NETWORK ARE CONSIDERED
  OBSERVED
  -----------------------------------------------------------------------------
*/
double
LSHKState(int K, struct group *part)
{
  struct group *g1=NULL, *g2=NULL;
  int n, nk;
  double H=0.0;
  int ng=0;      /* Number of non-empty groups */
  int nnod=0;  /* Number of nodes */
  int k;

  /* Go through all group pairs */
  g1 = part;
  while (g1->next != NULL) {
    g2 = g1;
    g1 = g1->next;
    if (g1->size > 0)
      ng++;
    nnod += g1->size;
    while ((g2 = g2->next) != NULL) {
      n = NG2GLinks(g1, g2);  /* The total number of observations */
      H += gsl_sf_lnfact(n + K - 1);
      for (k=0; k<K; k++) {
	nk = NWeightG2GLinks(g1, g2, (double)k);  /* The k ratings */
	H -= gsl_sf_lnfact(nk);
      }
    } /* end of loop over g2 */
  } /* end of loop over g1 */

  H -= gsl_sf_lnfact(nnod - ng);

  /* Done */
  return H;
}


/*
  -----------------------------------------------------------------------------
  Do a Monte Carlo step for the k-state link score. Fast version for
  sparse observation matrices.
  -----------------------------------------------------------------------------
*/
void
LSMCStepKState(int K,
	       int factor,
	       double *H,
	       struct node_gra **nlist,
	       struct group **glist,
	       struct group *part,
	       int nnod,
	       int *ng,
	       int ****N2G,
	       int ****G2G,
	       double *LogList, int LogListSize,
	       double *LogFactList, int LogFactListSize,
	       gsl_rng *gen)
{
  double dH;
  struct group *oldg, *newg, *g, *g2;
  int move;
  struct node_gra *node=NULL;
  struct node_lis *nei=NULL;
  int dice, n, nk;
  int oldgnum, newgnum, q;
  int i, ngroup=nnod;
  int j; 
  int move_in;
  int k;

  /* The moves */
  for (move=0; move<nnod*factor; move++) {

    /* The move */
    dice = floor(gsl_rng_uniform(gen) * (double)nnod);
    node = nlist[dice];
    oldgnum = node->inGroup;
    do {
      newgnum = floor(gsl_rng_uniform(gen) * (double)nnod);
    } while (newgnum == oldgnum);
    oldg = glist[oldgnum];
    newg = glist[newgnum];
    g = g2 = part;

    /* THE CHANGE OF ENERGY: OLD CONFIGURATION CONTRIBUTION */
    dH = 0.0;

    while ((g=g->next) != NULL) {
      /* if (g->size > 0) {  /\* group is not empty *\/ */
	/* old configuration, old group */
	n = 0;
	for (k=0; k<K; k++) {
	  n += (*G2G)[k][oldgnum][g->label];
	  nk = (*G2G)[k][oldgnum][g->label];
	  dH -= -FastLogFact(nk, LogFactList, LogFactListSize);
	}
	dH -= FastLogFact(n + K - 1, LogFactList, LogFactListSize);
	/* old configuration, new group */
	if (g->label != oldgnum) {
	  n = 0;
	  for (k=0; k<K; k++) {
	    n += (*G2G)[k][newgnum][g->label];
	    nk = (*G2G)[k][newgnum][g->label];
	    dH -= -FastLogFact(nk, LogFactList, LogFactListSize);
	  }
	  dH -= FastLogFact(n + K - 1, LogFactList, LogFactListSize);
	}
      /* } */
    }
    /* labeled-groups sampling correction */
    if ((oldg->size == 1) || (newg->size == 0)) {
      dH += FastLogFact(nnod - *ng, LogFactList, LogFactListSize);
    }

    /* Tentatively move the node to the new group and update G2G and N2G
       matrices */
    MoveNode(node, oldg, newg);
    for (i=0; i<ngroup; i++) {   /* update G2G links */
      for (k=0; k<K; k++) {
	(*G2G)[k][oldgnum][i] -= (*N2G)[k][node->num][i];
	(*G2G)[k][newgnum][i] += (*N2G)[k][node->num][i];
	(*G2G)[k][i][oldgnum] = (*G2G)[k][oldgnum][i];
	(*G2G)[k][i][newgnum] = (*G2G)[k][newgnum][i];
      }
    }
    nei = node->neig;            /* update N2G links */
    while ((nei = nei->next) != NULL) {
      for (k=0; k<K; k++) {
	if (nei->weight == (double)k) {
	  (*N2G)[k][nei->ref->num][oldgnum] -= 1;
	  (*N2G)[k][nei->ref->num][newgnum] += 1;
	}
      }
    }
    if (oldg->size == 0) /* update number of non-empty groups */
      (*ng) -= 1;
    if (newg->size == 1)
      (*ng) +=1;

    /* THE CHANGE OF ENERGY: NEW CONFIGURATION CONTRIBUTION */
    while ((g2=g2->next) != NULL) {
      /* if (g2->size > 0) {  /\* group is not empty *\/ */
	/* new configuration, old group */
	n = 0;
	for (k=0; k<K; k++) {
	  n += (*G2G)[k][oldgnum][g2->label];
	  nk = (*G2G)[k][oldgnum][g2->label];
	  dH += -FastLogFact(nk, LogFactList, LogFactListSize);
	}
	dH += FastLogFact(n + K - 1, LogFactList, LogFactListSize);
	/* new configuration, new group */
	if (g2->label != oldgnum) {
	  n = 0;
	  for (k=0; k<K; k++) {
	    n += (*G2G)[k][newgnum][g2->label];
	    nk = (*G2G)[k][newgnum][g2->label];
	    dH += -FastLogFact(nk, LogFactList, LogFactListSize);
	  }
	  dH += FastLogFact(n + K - 1, LogFactList, LogFactListSize);
	}
      /* } */
    }
    /* labeled-groups sampling correction */
    if ((oldg->size == 0) || (newg->size == 1)) {
      dH -= FastLogFact(nnod - *ng, LogFactList, LogFactListSize);
    }

    /* METROPOLIS ACCEPTANCE */
    if ((dH <= 0.0) || (gsl_rng_uniform(gen) < exp(-dH))) {
      /* accept move: update energy */
      *H += dH;
    }
    else {
      /* undo the move */
      MoveNode(node, newg, oldg);
      for (i=0; i<ngroup; i++) {   /* update G2G links */
	for (k=0; k<K; k++) {
	  (*G2G)[k][oldgnum][i] += (*N2G)[k][node->num][i];
	  (*G2G)[k][newgnum][i] -= (*N2G)[k][node->num][i];
	  (*G2G)[k][i][oldgnum] = (*G2G)[k][oldgnum][i];
	  (*G2G)[k][i][newgnum] = (*G2G)[k][newgnum][i];
	}
      }
      nei = node->neig;            /* update N2G links */
      while ((nei = nei->next) != NULL) {
	for (k=0; k<K; k++) {
	  if (nei->weight == (double)k) {
	    (*N2G)[k][nei->ref->num][oldgnum] += 1;
	    (*N2G)[k][nei->ref->num][newgnum] -= 1;
	  }
	}
      }
      if (oldg->size == 1) /* update number of non-empty groups */
	(*ng) += 1;
      if (newg->size == 0)
	(*ng) -= 1;
    }
  } /* Moves completed: done! */
  
  return;
}


/*
  -----------------------------------------------------------------------------
  Get the decorrelation step necessary to sample decorrelated
  partitions.
  -----------------------------------------------------------------------------
*/
int
LSGetDecorrelationStepKState(int K,
			     double *H,
			     struct node_gra **nlist,
			     struct group **glist,
			     struct group *part,
			     int nnod,
			     int *ng,
			     int ****N2G,
			     int ****G2G,
			     double *LogList, int LogListSize,
			     double *LogFactList, int LogFactListSize,
			     gsl_rng *gen,
			     char verbose_sw)
{
  struct group *partRef;
  int step, x1, x2;
  double y1, y2;
  double mutualInfo;
  int rep, nrep=10;
  double *decay, meanDecay, sigmaDecay, result;
  int norm=0;

  x2 = nnod / 5;
  if (x2 < 10)
    x2 = 10;
  x1 = x2 / 4;

  /* Get the nrep initial estimates */
  decay = allocate_d_vec(nrep);
  for (rep=0; rep<nrep; rep++) {
    switch (verbose_sw) {
    case 'q':
      break;
    default:
      fprintf(stderr, "#\n# Estimating decorrelation time (%d/%d)\n",
	      rep + 1, nrep);
      break;
    }
    partRef = CopyPartition(part);
    for (step=0; step<=x2; step++) {
      switch (verbose_sw) {
      case 'd':
	fprintf(stderr, "# %d / %d\n", step, x2);
	break;
      default:
	break;
      }
      LSMCStepKState(K, 1, H, nlist,
		     glist, part,
		     nnod, ng,
		     N2G,
		     G2G,
		     LogList, LogListSize,
		     LogFactList, LogFactListSize,
		     gen);
      if (step == x1) {
	y1 = MutualInformation(partRef, part, 0);
      }
    }
    y2 = MutualInformation(partRef, part, 0);
    decay[rep] = 2. * CalculateDecay(nnod, x1, y1, x2, y2);
    switch (verbose_sw) {
    case 'q':
      break;
    default:
      fprintf(stderr, "# Decorrelation time (estimate %d) = %g\n",
	      rep + 1, decay[rep]);
      break;
    }
    if (decay[rep] < 0.) {
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
    RemovePartition(partRef);
  }
  
  /* Get rid of bad estimates (Chauvenet criterion)  */
  meanDecay = mean(decay, nrep);
  sigmaDecay = stddev(decay, nrep);
  result = meanDecay * nrep;
  for (rep=0; rep<nrep; rep++) {
    if (fabs(decay[rep] - meanDecay) / sigmaDecay > 2) {
      result -= decay[rep];
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
  free_d_vec(decay);

  return (int)(result / norm + 0.5);
}


/*
  ---------------------------------------------------------------------
  
  ---------------------------------------------------------------------
*/
void
LSThermalizeMCKState(int K,
		     int decorStep,
		     double *H,
		     struct node_gra **nlist,
		     struct group **glist,
		     struct group *part,
		     int nnod,
		     int *ng,
		     int ****N2G,
		     int ****G2G,
		     double *LogList, int LogListSize,
		     double *LogFactList, int LogFactListSize,
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
      LSMCStepKState(K, decorStep, H, nlist,
		     glist, part,
		     nnod, ng,
		     N2G,
		     G2G,
		     LogList, LogListSize,
		     LogFactList, LogFactListSize,
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
  Return the score p(A_ij=k|A^O) of a collection {(i,j)} querySet of
  links, for a k-state system. The ratings are in a regular network
  with links (corresponding to observations) that have values 0 or
  1.
  -----------------------------------------------------------------------------
*/
double **
LSMultiLinkScoreKState(int K,
		       struct node_gra *observed,
		       int nIter,
		       gsl_rng *gen,
		       char verbose_sw,
		       int decorStep,
		       struct query ***querySet,
		       int *nquery)
{
  int nnod=CountNodes(observed);
  struct group *part=NULL;
  struct node_gra *p=NULL, *p2=NULL;
  struct node_gra **nlist=NULL;
  struct group **glist=NULL;
  struct group *lastg=NULL;
  double H;
  int iter;
  double **score;
  int i, j;
  int ***N2G;
  int ***G2G;
  int LogListSize = 5000;
  double *LogList=InitializeFastLog(LogListSize);
  int LogFactListSize = 10000;
  double *LogFactList=InitializeFastLogFact(LogFactListSize);
  int norm = 0;
  int dice;
  int n, nk;
  int q;
  int ng;
  FILE *outfile=NULL;
  int k, k2;
  int nq=0;

  /*
    PRELIMINARIES
  */
  N2G = (int ***) calloc(K, sizeof(int **));
  G2G = (int ***) calloc(K, sizeof(int **));

  /* Map nodes and groups to a list for faster access */
  fprintf(stderr, ">> Mapping nodes and groups to lists...\n");
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  glist = (struct group **) calloc(nnod, sizeof(struct group *));
  lastg = part = CreateHeaderGroup();
  p = observed;
  while ((p = p->next) != NULL) {
    nlist[p->num] = p;
    lastg = glist[p->num] = CreateGroup(lastg, p->num);
  }

  /* Create the queries: all pairs of nodes are queries */
  *nquery = (nnod * (nnod - 1)) / 2;
  (*querySet) = (struct query **) calloc(*nquery, sizeof(struct query *));
  p = observed;
  while ((p = p->next) != NULL) {
    p2 = p;
    while ((p2 = p2->next) != NULL) {
      (*querySet)[nq++] = CreateQuery(p, p2);
    }
  }

  /* Initialize scores */
  score = allocate_d_mat(K, *nquery);
  for (k=0; k<K; k++)
    for (q=0; q<*nquery; q++)
      score[k][q] = 0.0;

  /* Create initial partition (each node in a separate group) */
  fprintf(stderr, ">> Placing nodes in initial partition...\n");
  p = observed;
  ResetNetGroup(observed);
  while ((p = p->next) != NULL) {
    AddNodeToGroup(glist[p->num], p);
  }

  /* Get the initial group-to-group links matrix */
  fprintf(stderr, ">> Getting the initial group-to-group links matrix...\n");
  for (k=0; k<K; k++) {
    G2G[k] = allocate_i_mat(nnod, nnod);
    for (i=0; i<nnod; i++) {
      for (j=0; j<nnod; j++) {
	G2G[k][i][j] = G2G[k][j][i] = 
	  NWeightG2GLinks(glist[i], glist[j], (double)k);
      }
    }
  }

  /* Get the initial node-to-group links matrix */
  fprintf(stderr, ">> Getting the initial node-to-group links matrix...\n");
  for (k=0; k<K; k++) {
    N2G[k] = allocate_i_mat(nnod, nnod);
    for (i=0; i<nnod; i++) {
      for (j=0; j<nnod; j++) {
	N2G[k][i][j] = NWeightLinksToGroup(nlist[i], glist[j], (double)k);
      }
    }
  }

  /* Get the initial number of non-empty groups */
  ng = NNonEmptyGroups(part);

  /*
    GET READY FOR THE SAMPLING
  */
  H = LSHKState(K, part);

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
    decorStep = LSGetDecorrelationStepKState(K, &H,
  					     nlist,
  					     glist,
  					     part,
  					     nnod,
  					     &ng,
  					     &N2G,
  					     &G2G,
  					     LogList, LogListSize,
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
  LSThermalizeMCKState(K, decorStep, &H,
  		       nlist, glist,
  		       part,
  		       nnod, &ng,
  		       &N2G,
  		       &G2G,
  		       LogList, LogListSize,
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
    LSMCStepKState(K, decorStep, &H, nlist,
		   glist, part,
		   nnod, &ng,
		   &N2G,
		   &G2G,
		   LogList, LogListSize,
		   LogFactList, LogFactListSize,
		   gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf\n", iter, H, LSHKState(K, part));
      /* FPrintPartition(stderr, part, 0); */
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
      for (k=0; k<K; k++) {
	for (q=0; q<*nquery; q++) {
	  score[k][q] = 0.0;
	}
      }
    }

    /* Update normalization */
    norm += 1;

    /* Update the scores */
    for (k=0; k<K; k++) {
      for (q=0; q<*nquery; q++) {
	nk = G2G[k][(*querySet)[q]->n1->inGroup][(*querySet)[q]->n2->inGroup];
	n = 0;
	for (k2=0; k2<K; k2++)
	  n += G2G[k2][(*querySet)[q]->n1->inGroup][(*querySet)[q]->n2->inGroup];
	score[k][q] += (float)(nk + 1) / (float)(n + K);
      }
    }

    /* Output temporary scores and partitions */
    if (iter % 100 == 0) {
      switch (verbose_sw) {
      case 'q':
	break;
      default:
 	outfile = fopen("scores.tmp", "w");
	for (q=0; q<*nquery; q++) {
	  fprintf(outfile, "%s %s",
		  (((*querySet)[q])->n1)->label,
		  (((*querySet)[q])->n2)->label);
	  for (k=0; k<K; k++) {
	    fprintf(outfile, " %lf",
		    score[k][q]/(double)norm);
	  }
	  fprintf(outfile, "\n");
	}
	fclose(outfile);
	outfile = fopen("part.tmp", "w");
	FPrintPartition(outfile, part, 1);
	fclose(outfile);
	break;
      }
    }
  }  /* End of iter loop */

  /* Normalize the scores */
  for (k=0; k<K; k++)
    for (q=0; q<*nquery; q++)
      score[k][q] /= (double)norm;
  
  /* Free dynamically allocated memory */
  RemovePartition(part);
  free(glist);
  free(nlist);
  for (k=0; k<K; k++) {
    free_i_mat(G2G[k], nnod);
    free_i_mat(N2G[k], nnod);
  }
  free(G2G);
  free(N2G);
  FreeFastLog(LogList);
  FreeFastLogFact(LogFactList);

  /* Done */
  return score;
}





/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Gibbs sampling
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Do a Gibbs Monte Carlo step for the prediction of missing links.
  ---------------------------------------------------------------------
*/
void
GibbsLinkScoreStep(double *H,
		   double linC,
		   struct node_gra **nlist,
		   struct group **glist,
		   struct group *part,
		   int nnod,
		   int **G2G,
		   int *n2gList,
		   double **LogChooseList,
		   int LogChooseListSize,
		   double *LogFactList, int LogFactListSize,
		   gsl_rng *gen)
{
  double *dH;
  struct group *g, *oldg, *newg;
  int oldgnum, newgnum, n2g, oldg2g, newg2g, ng, noldg, nnewg;
  struct node_gra *node;
  int n2oldg, n2newg, newg2oldg, oldg2oldg, newg2newg;
  int r, l;
  int i, j, ngroup=nnod;
  int move;
  int effectng;
  /* Shortlist of groups. Groups not in this list are empty for
     sure. Note that, during a step, we may have to add groups to this
     list (groups that were initially empty but got filled) but we
     never remove groups from it (even if they loose all their nodes
     during a step). slgXsize is the size of slgX. */
  int *slg=NULL;
  int slgn, slgsize, isinlist;
  int newgn, target;
  double dice, norm, dHempty, cum;
  
  /* Preliminaries */
  dH = allocate_d_vec(nnod);
  slg = allocate_i_vec(nnod);  // Shortlist of groups
  slgsize = 0;
  g = part;
  while ((g=g->next) != NULL)
    if (g->size > 0)
      slg[slgsize++] = g->label;

  /* Steps */
  for (move=0; move<nnod; move++) {
    node = nlist[move];
    oldgnum = node->inGroup;
    oldg = glist[oldgnum];
    norm = 0.0;

    /* Loop over destination groups */
    for (newgn=0; newgn<slgsize; newgn++) {
      newgnum = slg[newgn];
      newg = glist[newgnum];

      if (newgnum != oldgnum) {	
	/* Calculate the change of energy */
	dH[newgnum] = 0.0;
	noldg = oldg->size;
	nnewg = newg->size;
	if (noldg == 1)  /* number of groups would decrease by one */ 
	  dH[newgnum] -= linC;
	if (nnewg == 0)  /* number of groups would increase by one */
	  dH[newgnum] += linC;
	n2oldg = NLinksToGroup(node, oldg);
	n2newg = NLinksToGroup(node, newg);
	newg2oldg = NG2GLinks(newg, oldg);
	oldg2oldg = oldg->inlinks;
	newg2newg = newg->inlinks;
	effectng = 0;
	for (slgn=0; slgn<slgsize; slgn++) {
	  g = glist[slg[slgn]];
	  if (g->size > 0) {  /* group is not empty */
	    effectng++;
	    n2gList[g->label] = NLinksToGroup(node, g);
	    if (g->label == oldg->label) {
	      /* old conf, oldg-oldg */
	      r = noldg * (noldg - 1) / 2;
	      l = oldg2oldg;
	      dH[newgnum] -= log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /* new conf, oldg-olg */
	      r = (noldg - 1) * (noldg - 2) / 2;
	      l = oldg2oldg - n2oldg;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /* old conf, newg-oldg */
	      r = nnewg * noldg;
	      l = newg2oldg;
	      dH[newgnum] -= log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /* new conf, newg-oldg */
	      r = (nnewg + 1) * (noldg - 1);
	      l = newg2oldg + n2oldg - n2newg;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	    }
	    else if (g->label == newg->label) {
	      /* old conf, newg-newg */
	      r = nnewg * (nnewg - 1) / 2;
	      l = newg2newg;
	      dH[newgnum] -= log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /* new conf, newg-olg */
	      r = (nnewg + 1) * nnewg / 2;
	      l = newg2newg + n2newg;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	    }
	    else {
	      n2g = n2gList[g->label];
	      oldg2g = G2G[oldg->label][g->label];
	      newg2g = G2G[newg->label][g->label];
	      ng = g->size;
	      /* old conf, oldg-g */
	      r = noldg * ng;
	      l = oldg2g;
	      dH[newgnum] -= log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /* new conf, oldg-g */
	      r = (noldg - 1) * ng;
	      l = oldg2g - n2g;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /* old conf, newg-g */
	      r = nnewg * ng;
	      l = newg2g;
	      dH[newgnum] -= log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /* new conf, newg-g */
	      r = (nnewg + 1) * ng;
	      l = newg2g + n2g;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	    }
	  }
	  else { /* group is empty */
	    n2gList[g->label] = 0.0;
	  }
	}

	/* Labeled-groups sampling correction */
	if (noldg == 1 && nnewg != 0) { // effectng would decrease by one
	  dH[newgnum] -= -FastLogFact(ngroup - effectng,
				      LogFactList, LogFactListSize);
	  dH[newgnum] += -FastLogFact(ngroup - (effectng-1),
				      LogFactList, LogFactListSize);
	}
	else if (noldg != 1 && nnewg == 0) {  // effectng would increase by one
	  dH[newgnum] -= -FastLogFact(ngroup - effectng,
				      LogFactList, LogFactListSize);
	  dH[newgnum] += -FastLogFact(ngroup - (effectng+1),
				      LogFactList, LogFactListSize);
	}
      }

      else { // oldg and newg are the same: nothing changes
	dH[newgnum] = 0.0;
      }

      norm += exp(-dH[newgnum]);
    }

    /** Calculate the change of energy to go to an empty group **/
    dHempty = 0.0;
    nnewg = 0;
    if (noldg == 1)  /* number of groups would decrease by one */ 
      dHempty -= linC;
    if (nnewg == 0)  /* number of groups would increase by one */
      dHempty += linC;
    n2newg = 0;
    newg2oldg = 0;
    newg2newg = 0;

    for (slgn=0; slgn<slgsize; slgn++) {
      g = glist[slg[slgn]];
      if (g->size > 0) {  /* group is not empty */
	if (g->label == oldg->label) {
	  /* old conf, oldg-oldg */
	  r = noldg * (noldg - 1) / 2;
	  l = oldg2oldg;
	  dHempty -= log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /* new conf, oldg-olg */
	  r = (noldg - 1) * (noldg - 2) / 2;
	  l = oldg2oldg - n2oldg;
	  dHempty += log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /* old conf, newg-oldg */
	  r = nnewg * noldg;
	  l = newg2oldg;
	  dHempty -= log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /* new conf, newg-oldg */
	  r = (nnewg + 1) * (noldg - 1);
	  l = newg2oldg + n2oldg - n2newg;
	  dHempty += log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	}
	else {
	  n2g = n2gList[g->label];
	  oldg2g = G2G[oldg->label][g->label];
	  newg2g = 0;
	  ng = g->size;
	  /* old conf, oldg-g */
	  r = noldg * ng;
	  l = oldg2g;
	  dHempty -= log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /* new conf, oldg-g */
	  r = (noldg - 1) * ng;
	  l = oldg2g - n2g;
	  dHempty += log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /* old conf, newg-g */
	  r = nnewg * ng;
	  l = newg2g;
	  dHempty -= log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /* new conf, newg-g */
	  r = (nnewg + 1) * ng;
	  l = newg2g + n2g;
	  dHempty += log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	}
      }
    }

    /* Labeled-groups sampling correction */
    if (noldg == 1 && nnewg != 0) { // effectng would decrease by one
      dHempty -= -FastLogFact(ngroup - effectng,
			      LogFactList, LogFactListSize);
      dHempty += -FastLogFact(ngroup - (effectng-1),
			      LogFactList, LogFactListSize);
    }
    else if (noldg != 1 && nnewg == 0) {  // effectng would increase by one
      dHempty -= -FastLogFact(ngroup - effectng,
			      LogFactList, LogFactListSize);
      dHempty += -FastLogFact(ngroup - (effectng+1),
			      LogFactList, LogFactListSize);
    }

    norm += exp(-dHempty) * (double)(nnod - slgsize);
    
    /** CHOOSE THE MOVE **/
    dice = norm * gsl_rng_uniform(gen);

    if (dice > (norm - (double)(nnod - slgsize) * exp(-dHempty))) {
      /* Select an empty group */
      newg = GetEmptyGroup(part);
      newgnum = newg->label;
      dH[newgnum] = dHempty;
    }
    else {
      /* Select group in shortlist */
      target = 0;
      cum = 0.0;
      while (cum < dice) {
	cum += exp(-dH[slg[target++]]);
      }
      newgnum = slg[target - 1];
      newg = glist[newgnum];
    }
    
    /* MAKE THE MOVE AND UPDATE MATRICES */
    MoveNode(node, oldg, newg);
    *H += dH[newgnum];
    
    /* update G2G */
    for (slgn=0; slgn<slgsize; slgn++) {
      G2G[slg[slgn]][oldgnum] -= n2gList[slg[slgn]];
      G2G[oldgnum][slg[slgn]] = G2G[slg[slgn]][oldgnum];
      G2G[slg[slgn]][newgnum] += n2gList[slg[slgn]];
      G2G[newgnum][slg[slgn]] = G2G[slg[slgn]][newgnum];
    }
    
    /* Add newg to shortlist (if it's not there already!) */
    isinlist = 0;
    for (slgn=0; slgn<slgsize; slgn++)
      if (newgnum == slg[slgn])
	isinlist = 1;
    if (isinlist == 0)
      slg[slgsize++] = newgnum;
  }  /* nnod moves completed: done! */

  free_i_vec(slg);
  free_d_vec(dH);
}


/*
  ---------------------------------------------------------------------
  
  ---------------------------------------------------------------------
*/
void
GibbsThermalizeLinkScore(double *H,
			 double linC,
			 struct node_gra **nlist,
			 struct group **glist,
			 struct group *part,
			 int nnod,
			 int **G2G,
			 int *n2gList,
			 double **LogChooseList,
			 int LogChooseListSize,
			 double *LogFactList, int LogFactListSize,
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
      GibbsLinkScoreStep(H, linC, nlist, glist, part,
			 nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
			 LogFactList, LogFactListSize,
			 gen);
      switch (verbose_sw) {
      case 'q':
	break;
      case 'd':
	fprintf(stderr, "%lf %lf %d\n", *H, PartitionH(part, linC),
		NNonEmptyGroups(part));
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
    if ((HMean0 - HStd0 / sqrt(nrep)) - (HMean1 + HStd1 / sqrt(nrep))
	< EPSILON) {
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
  ---------------------------------------------------------------------
  Predict the links missing in a network. The algorithm returns a
  matrix of scores for all links.
  ---------------------------------------------------------------------
*/
double **
GibbsLinkScore(struct node_gra *net,
	       double linC,
	       int nIter,
	       gsl_rng *gen,
	       char verbose_sw)
{
  int nnod=CountNodes(net);
  struct group *part=NULL;
  struct node_gra *p=NULL, *node=NULL;
  struct node_gra **nlist=NULL;
  struct group **glist=NULL;
  struct group *lastg=NULL;
  double H;
  int iter, decorStep;
  double **predA=NULL;
  int i, j;
  int **G2G=NULL;
  int *n2gList=NULL;
  int LogChooseListSize = 10000;
  double **LogChooseList=InitializeFastLogChoose(LogChooseListSize);
  struct node_lis *p1=NULL, *p2=NULL;
  double contrib;
  int dice;
  int norm = 0;
  int r, l;
  double mutualInfo;
  int LogFactListSize = 10000;
  double *LogFactList=InitializeFastLogFact(LogFactListSize);

  /*
    PRELIMINARIES
  */
  /* Initialize the predicted adjacency matrix */
  predA = allocate_d_mat(nnod, nnod);
  for (i=0; i<nnod; i++)
    for (j=0; j<nnod; j++)
      predA[i][j] = 0.0;

  /* Map nodes and groups to a list for faster access */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  glist = (struct group **) calloc(nnod, sizeof(struct group *));
  lastg = part = CreateHeaderGroup();
  p = net;
  while ((p = p->next) != NULL) {
    nlist[p->num] = p;
    lastg = glist[p->num] = CreateGroup(lastg, p->num);
  }

  /* Place nodes in random partitions */
  p = net;
  ResetNetGroup(net);
  while ((p = p->next) != NULL) {
    dice = floor(gsl_rng_uniform(gen) * (double)nnod);
    /* AddNodeToGroup(glist[dice], p); */
    AddNodeToGroup(glist[p->num], p);
  }

  /* Get the initial group-to-group links matrix */
  G2G = allocate_i_mat(nnod, nnod);
  n2gList = allocate_i_vec(nnod);
  for (i=0; i<nnod; i++) {
    G2G[i][i] = glist[i]->inlinks;
    for (j=i+1; j<nnod; j++) {
      G2G[i][j] = G2G[j][i] = NG2GLinks(glist[i], glist[j]);
    }
  }

  /* Initial energy */
  H = PartitionH(part, linC);

  /*
    GET READY FOR THE SAMPLING
  */

  /* Thermalization */
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "#\n#\n# THERMALIZING\n");
    fprintf(stderr, "# ------------\n");
    break;
  }
  GibbsThermalizeLinkScore(&H, linC, nlist, glist, part,
			   nnod, G2G, n2gList,
			   LogChooseList, LogChooseListSize,
			   LogFactList, LogFactListSize,
			   gen, verbose_sw);
  
  /*
    SAMPLIN' ALONG
  */
  /* Unless we are in debug mode, reset the origin of energies */
  switch (verbose_sw) {
  case 'd':
    break;
  default:
    /* H = 0; */
    break;
  }

  /* Do the MC Steps */
  for (iter=0; iter<nIter; iter++) {
    GibbsLinkScoreStep(&H, linC, nlist, glist, part,
		       nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
		       LogFactList, LogFactListSize,
		       gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf %d\n", iter, H, PartitionH(part, linC),
	      NNonEmptyGroups(part));
      /* FPrintPartition(stderr, part, 0); */
      break;
    }

    /* Update partition function */
    norm += 1;

    /* Update the predicted adjacency matrix by going through all
       group pairs */
    for (i=0; i<nnod; i++) {
      if (glist[i]->size > 0) {
	
	/* update the within-group pairs */
	r = glist[i]->size * (glist[i]->size - 1) / 2;
	l = glist[i]->inlinks;
	contrib = (float)(l + 1) / (float)(r + 2);
	p1 = glist[i]->nodeList;
	while ((p1 = p1->next) != NULL) {
	  p2 = p1;
	  while ((p2 = p2->next) != NULL) {
	    predA[p1->node][p2->node] += contrib;
	    predA[p2->node][p1->node] += contrib;
	  }
	}
      
	/* update the between-group pairs */
	for (j=i+1; j<nnod; j++) {
	  if (glist[j]->size > 0) {
	    l = G2G[i][j];
	    r = glist[i]->size * glist[j]->size;
	    contrib = (float)(l + 1) / (float)(r + 2);
	    p1 = glist[i]->nodeList;
	    while ((p1 = p1->next) != NULL) {
	      p2 = glist[j]->nodeList;
	      while ((p2 = p2->next) != NULL) {
		predA[p1->node][p2->node] += contrib;
		predA[p2->node][p1->node] += contrib;
	      }
	    }
	  }
	}
      }
    } /* Done updating adjacency matrix */

  }  /* End of iter loop */

  /* Normalize the predicted adjacency matrix */
  for (i=0; i<nnod; i++) {
    for (j=0; j<nnod; j++) {
      predA[i][j] /= (double)norm;
    }
  }

  /* Done */
  RemovePartition(part);
  free(glist);
  free(nlist);
  free_i_mat(G2G, nnod);
  free_i_vec(n2gList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  FreeFastLogFact(LogFactList);

  return predA;
}
