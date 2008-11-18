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

#include "prng.h"

#include "tools.h"
#include "graph.h"
#include "modules.h"

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
		 gsl_vector * f)
{
  const double a = gsl_vector_get(params, 0);
  const double b = gsl_vector_get(params, 1);

  double x0 = ((struct pair *)points)->x0;
  double y0 = ((struct pair *)points)->y0;
  double x1 = ((struct pair *)points)->x1;
  double y1 = ((struct pair *)points)->y1;
  
  const double r0 = y0 - a - (1. - a) * exp(-x0 / b);
  const double r1 = y1 - a - (1. - a) * exp(-x1 / b);
  
  gsl_vector_set (f, 0, r0);
  gsl_vector_set (f, 1, r1);
  
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
  
  fprintf(stderr, "# GSL status = %s\n", gsl_strerror(status));
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

  H += linC * NNonEmptyGroups(part);

  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
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

  return H;
}

/*
  ---------------------------------------------------------------------
  Do a Monte Carlo step for the prediction of missing links. Each
  "step" involves nnod independent attempts to move a single node.
  ---------------------------------------------------------------------
*/
void
MissingLinksMCStep(int factor,
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
		   struct prng *gen)
{
  double dH;
  struct group *g, *oldg, *newg;
  int dice, oldgnum, newgnum, n2g, oldg2g, newg2g, ng, noldg, nnewg;
  struct node_gra *node;
  int n2oldg, n2newg, newg2oldg, oldg2oldg, newg2newg;
  int r, l;
  int i, j;
  int move;

  for (move=0; move<nnod*factor; move++) {

    /* Choose node and destination group */
    dice = floor(prng_get_next(gen) * (double)nnod);
    node = nlist[dice];
    oldgnum = node->inGroup;
    do {
      newgnum = floor(prng_get_next(gen) * (double)nnod);
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
    g = part;
    while ((g = g->next) != NULL) {
      if (g->size > 0) {  /* group is not empty */
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

    /* Metropolis rule */
    if ((dH <= 0.0) || (prng_get_next(gen) < exp(-dH))) {
      
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
    }
  } /* nnod moves completed: done! */
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
		     struct prng *gen)
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
    fprintf(stderr, "#\n# Estimating decorrelation time (%d/%d)\n",
	    rep + 1, nrep);
    partRef = CopyPartition(part);
    for (step=0; step<=x2; step++) {
      MissingLinksMCStep(1, H, linC, nlist, glist, part,
			 nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
			 gen);
      if (step == x1)
	y1 = MutualInformation(partRef, part);
    }
    y2 = MutualInformation(partRef, part);
    RemovePartition(partRef);
    decay[rep] = 2. * CalculateDecay(nnod, x1, y1, x2, y2);
    fprintf(stderr, "# Decorrelation time (estimate %d) = %g\n",
	    rep + 1, decay[rep]);
    if (decay[rep] < 0) {
      fprintf(stderr, "#\tignoring...\n");
      rep--;
    }
  }
  fprintf(stderr, "#\n");
  
  /* Get rid of bad estimates (Chauvenet criterion)  */
  meanDecay = mean(decay, nrep);
  sigmaDecay = stddev(decay, nrep);
  result = meanDecay * nrep;
  for (rep=0; rep<nrep; rep++) {
    if (fabs(decay[rep] - meanDecay) / sigmaDecay > 2) {
      result -= decay[rep];
      fprintf(stderr, "# Disregarding estimate %d\n", rep + 1);
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
ThermalizeMissingLinkMC(int decorStep,
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
			struct prng *gen)
{
  double HMean0=1.e10, HStd0=1.e-10, HMean1, HStd1, *Hvalues;
  int rep, nrep=20;
  int equilibrated=0;

  Hvalues = allocate_d_vec(nrep);

  do {
    
    /* MC steps */
    for (rep=0; rep<nrep; rep++) {
      MissingLinksMCStep(decorStep, H, linC, nlist, glist, part,
			 nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
			 gen);
      fprintf(stderr, "%lf\n", *H);
      Hvalues[rep] = *H;
    }

    /* Check for equilibration */
    HMean1 = mean(Hvalues, nrep);
    HStd1 = stddev(Hvalues, nrep);
    if (HMean0 - HStd0 / sqrt(nrep) < HMean1 + HStd1 / sqrt(nrep)) {
      fprintf(stderr, "#\tequilibrated (%d/5) H=%lf\n",
	      ++equilibrated, HMean1);
    }
    else {
      fprintf(stderr, "#\tnot equilibrated yet H0=%g+-%g H1=%g+-%g\n",
	      HMean0, HStd0 / sqrt(nrep), HMean1, HStd1 / sqrt(nrep));
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
MissingLinks(struct node_gra *net, double linC, int nIter, struct prng *gen)
{
  int nnod=CountNodes(net);
  struct group *part;
  struct node_gra *p, *node;
  struct node_gra **nlist;
  struct group **glist;
  struct group *lastg;
  double H;
  int iter, decorStep;
  double **predA, Z=0.0;
  int i, j;
  int **G2G;
  int *n2gList;
  int LogChooseListSize = 500;
  double **LogChooseList=InitializeFastLogChoose(LogChooseListSize);
  struct node_lis *p1, *p2;
  double weight, contrib;
  int dice;
  int r, l;
  double mutualInfo;

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
  while ((p = p->next) != NULL) {
    dice = floor(prng_get_next(gen) * (double)nnod);
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
  fprintf(stderr, "# CALCULATING DECORRELATION TIME\n");
  fprintf(stderr, "# ------------------------------\n");
  decorStep = GetDecorrelationStep(&H, linC, nlist, glist, part,
				   nnod, G2G, n2gList,
				   LogChooseList, LogChooseListSize,
				   gen);

  /* Thermalization */
  fprintf(stderr, "#\n#\n# THERMALIZING\n");
  fprintf(stderr, "# ------------\n");
  ThermalizeMissingLinkMC(decorStep, &H, linC, nlist, glist, part,
			  nnod, G2G, n2gList,
			  LogChooseList, LogChooseListSize,
			  gen);
  
  /*
    SAMPLIN' ALONG
  */
  H = 0; /* Reset the origin of energies to avoid huge exponentials */
  for (iter=0; iter<nIter; iter++) {
    MissingLinksMCStep(decorStep, &H, linC, nlist, glist, part,
		       nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
		       gen);
    fprintf(stderr, "%d %lf\n", iter, H);
/*     fprintf(stderr, "%d %lf %lf\n", iter, H, PartitionH(part, linC)); */

    /* Update partition function */
    weight = exp(-H);
    Z += weight;

    /* Update the predicted adjacency matrix by going through all
       group pairs */
    for (i=0; i<nnod; i++) {
      if (glist[i]->size > 0) {
	
	/* update the within-group pairs */
	r = glist[i]->size * (glist[i]->size - 1) / 2;
	l = glist[i]->inlinks;
	contrib = weight * (float)(l + 1) / (float)(r + 2);
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
	    contrib = weight * (float)(l + 1) / (float)(r + 2);
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
      predA[i][j] /= Z;
    }
  }

  /* Done */
  free(glist);
  free(nlist);
  free_i_mat(G2G, nnod);
  free_i_vec(n2gList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  
  return predA;
}

/*
  ---------------------------------------------------------------------
  Return the network reliability score
  ---------------------------------------------------------------------
*/
double
NetworkReliability(struct node_gra *net, struct prng *gen)
{
  struct node_gra **nlist;
  int i, j, nnod=CounNodes(net);
  struct node_gra *p;
  double score=0.0;
  double **pairScore;

  /* Get the link score */
  pairScore = MissingLinks(net, 0.0, 10000, gen);

  /* Map nodes to a list for faster access */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  p = net;
  while ((p = p->next) != NULL) {
    nlist[p->num] = p;
  }

  /* Calculate the reliability score */
  for (i=0; i<nnod; i++) {
    for (j=0; j<nnod; j++) {
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
  return score;
}
