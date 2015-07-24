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
  Sparse networks (beta priors on q)
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
/*
  ---------------------------------------------------------------------
  Partition H
  ---------------------------------------------------------------------
*/
double
SparsePartitionH(struct group *part, double a, double b)
{
  struct group *g1=part, *g2;
  double r, l, H=0.0;
  int nnod=0, ng=0;

  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
      ng++;
    }
    nnod += g1->size;
    r = g1->size * (g1->size - 1) / 2;
    l = g1->inlinks;
    H += gsl_sf_lngamma(r + a + b) - gsl_sf_lngamma(l + a) -	\
      gsl_sf_lngamma(r - l + b);
    g2 = g1;
    while ((g2 = g2->next) != NULL) {
      r = g1->size * g2->size;
      l = NG2GLinks(g1, g2);
      H += gsl_sf_lngamma(r + a + b) - gsl_sf_lngamma(l + a) -	\
	gsl_sf_lngamma(r - l + b);
    }
  }

  H -= gsl_sf_lnfact(nnod - ng);
  
  return H;
}



/*
  ---------------------------------------------------------------------
  Do a Gibbs Monte Carlo step for the prediction of missing links.
  ---------------------------------------------------------------------
*/
void
SparseGibbsLinkScoreStep(double *H,
			 struct node_gra **nlist,
			 struct group **glist,
			 struct group *part,
			 int nnod,
			 int **G2G,
			 int *n2gList,
			 double *LogGammaListA,
			 double *LogGammaListB,
			 double *LogGammaListAB,
			 int LogGammaListSize,
			 double *LogFactList, int LogFactListSize,
			 double a, double b,
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
	  }
	  n2gList[g->label] = NLinksToGroup(node, g);
	  if (g->label == oldg->label) {
	    /* old conf, oldg-oldg */
	    r = noldg * (noldg - 1) / 2;
	    l = oldg2oldg;
	    /* dH[newgnum] -= log(r + 1) + LogChoose(r, l); */
	    dH[newgnum] -= FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	    /* new conf, oldg-olg */
	    r = (noldg - 1) * (noldg - 2) / 2;
	    l = oldg2oldg - n2oldg;
	    dH[newgnum] += FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	    /* old conf, newg-oldg */
	    r = nnewg * noldg;
	    l = newg2oldg;
	    dH[newgnum] -= FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	    /* new conf, newg-oldg */
	    r = (nnewg + 1) * (noldg - 1);
	    l = newg2oldg + n2oldg - n2newg;
	    dH[newgnum] += FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	  }
	  else if (g->label == newg->label) {
	    /* old conf, newg-newg */
	    r = nnewg * (nnewg - 1) / 2;
	    l = newg2newg;
	    dH[newgnum] -= FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	    /* new conf, newg-olg */
	    r = (nnewg + 1) * nnewg / 2;
	    l = newg2newg + n2newg;
	    dH[newgnum] += FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	  }
	  else {
	    n2g = n2gList[g->label];
	    oldg2g = G2G[oldg->label][g->label];
	    newg2g = G2G[newg->label][g->label];
	    ng = g->size;
	    /* old conf, oldg-g */
	    r = noldg * ng;
	    l = oldg2g;
	    dH[newgnum] -= FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	    /* new conf, oldg-g */
	    r = (noldg - 1) * ng;
	    l = oldg2g - n2g;
	    dH[newgnum] += FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	    /* old conf, newg-g */
	    r = nnewg * ng;
	    l = newg2g;
	    dH[newgnum] -= FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	    /* new conf, newg-g */
	    r = (nnewg + 1) * ng;
	    l = newg2g + n2g;
	    dH[newgnum] += FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
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
    n2newg = 0;
    newg2oldg = 0;
    newg2newg = 0;

    for (slgn=0; slgn<slgsize; slgn++) {
      g = glist[slg[slgn]];
      if (g->label == oldg->label) {
	/* old conf, oldg-oldg */
	r = noldg * (noldg - 1) / 2;
	l = oldg2oldg;
	dHempty -= FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	/* new conf, oldg-olg */
	r = (noldg - 1) * (noldg - 2) / 2;
	l = oldg2oldg - n2oldg;
	dHempty += FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	/* old conf, newg-oldg */
	r = nnewg * noldg;
	l = newg2oldg;
	dHempty -= FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	/* new conf, newg-oldg */
	r = (nnewg + 1) * (noldg - 1);
	l = newg2oldg + n2oldg - n2newg;
	dHempty += FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
      }
      else {
	n2g = n2gList[g->label];
	oldg2g = G2G[oldg->label][g->label];
	newg2g = 0;
	ng = g->size;
	/* old conf, oldg-g */
	r = noldg * ng;
	l = oldg2g;
	dHempty -= FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	/* new conf, oldg-g */
	r = (noldg - 1) * ng;
	l = oldg2g - n2g;
	dHempty += FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	/* old conf, newg-g */
	r = nnewg * ng;
	l = newg2g;
	dHempty -= FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
	/* new conf, newg-g */
	r = (nnewg + 1) * ng;
	l = newg2g + n2g;
	dHempty += FastLogGamma(r, LogGammaListAB, LogGammaListSize, a+b) - FastLogGamma(l, LogGammaListA, LogGammaListSize, a) - FastLogGamma(r-l, LogGammaListB, LogGammaListSize, b);
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
      /* fprintf(stderr, "Move not in shortlist\n"); */
      /* Select an empty group */
      newg = GetEmptyGroup(part);
      newgnum = newg->label;
      dH[newgnum] = dHempty;
    }
    else {
      /* fprintf(stderr, "Move in shortlist\n"); */
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
SparseGibbsThermalizeLinkScore(double *H,
			       struct node_gra **nlist,
			       struct group **glist,
			       struct group *part,
			       int nnod,
			       int **G2G,
			       int *n2gList,
			       double *LogGammaListA,
			       double *LogGammaListB,
			       double *LogGammaListAB,
			       int LogGammaListSize,
			       double *LogFactList, int LogFactListSize,
			       double betaA, double betaB,
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
      SparseGibbsLinkScoreStep(H, nlist, glist, part,
			       nnod, G2G, n2gList,
			       LogGammaListA, LogGammaListB, LogGammaListAB,
			       LogGammaListSize,
			       LogFactList, LogFactListSize,
			       betaA, betaB,
			       gen);
      switch (verbose_sw) {
      case 'q':
	break;
      case 'd':
	fprintf(stderr, "%lf %lf %d\n",
		*H,
		SparsePartitionH(part, betaA, betaB),
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
SparseGibbsLinkScore(struct node_gra *net,
		     int nIter,
		     gsl_rng *gen,
		     char verbose_sw)
{
  int nnod=CountNodes(net), nlink=TotalNLinks(net, 1);
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
  struct node_lis *p1=NULL, *p2=NULL;
  double contrib;
  int dice;
  int norm = 0;
  int r, l;
  double mutualInfo;
  int LogGammaListSize = 10000;
  double *LogGammaListA=NULL, *LogGammaListB=NULL, *LogGammaListAB=NULL;
  int LogFactListSize = 10000;
  double *LogFactList=NULL;
  double betaA, betaB;
  double density;
  double Href;

  /*
    PRELIMINARIES
  */
  /* Initialize the parameters of the q priors */
  density = (2.0 * nlink / (double)(nnod * (nnod - 1)));
  betaB = .5;
  betaA = density * betaB / (1.0 - density);

  /* Initialize the fast Log(Gamma) arrays */
  LogGammaListA = InitializeFastLogGamma(LogGammaListSize, betaA);
  LogGammaListB = InitializeFastLogGamma(LogGammaListSize, betaB);
  LogGammaListAB = InitializeFastLogGamma(LogGammaListSize, betaA + betaB);
  LogFactList = InitializeFastLogFact(LogFactListSize);

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
  H = SparsePartitionH(part, betaA, betaB);

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
  SparseGibbsThermalizeLinkScore(&H, nlist, glist, part,
				 nnod, G2G, n2gList,
				 LogGammaListA, LogGammaListB, LogGammaListAB,
				 LogGammaListSize,
				 LogFactList, LogFactListSize,
				 betaA, betaB,
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
    SparseGibbsLinkScoreStep(&H, nlist, glist, part,
			     nnod, G2G, n2gList,
			     LogGammaListA, LogGammaListB, LogGammaListAB,
			     LogGammaListSize,
			     LogFactList, LogFactListSize,
			     betaA, betaB, gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf %d\n", iter, H, SparsePartitionH(part, betaA, betaB),
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
	contrib = (float)(l + betaA) / (float)(r + betaA + betaB);
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
	    contrib = (float)(l + betaA) / (float)(r + betaA + betaB);
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

  Href = SparsePartitionH(part, betaA, betaB);
  fprintf(stderr, "# Final H check: %g %g %g\n", H, Href, H-Href);
    
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
  FreeFastLogGamma(LogGammaListA);
  FreeFastLogGamma(LogGammaListB);
  FreeFastLogGamma(LogGammaListAB);
  FreeFastLogFact(LogFactList);

  return predA;
}
