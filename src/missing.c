/*
  missing.c
  $LastChangedDate: 2008-10-22 17:39:18 -0500 (Wed, 22 Oct 2008) $
  $Revision: 134 $
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include "prng.h"

#include "tools.h"
#include "graph.h"
#include "modules.h"

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Missing links
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
/*
  ---------------------------------------------------------------------
  Partition H
  ---------------------------------------------------------------------
*/
double
PartitionH(struct group *part)
{
  struct group *g1=part, *g2;
  double r, l, H=0.0;

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
  Do a Monte Carlo step for the prediction of missing links
  ---------------------------------------------------------------------
*/
void
MissingLinksMCStep(double *H,
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
}


/*
  ---------------------------------------------------------------------
  Predict the links missing in a network. The algorithm returns a
  matrix of scores for all links.
  ---------------------------------------------------------------------
*/
double **
MissingLinks(struct node_gra *net, struct prng *gen)
{
  int nnod=CountNodes(net);
  struct group *part, *partRef;
  struct node_gra *p, *node;
  struct node_gra **nlist;
  struct group **glist;
  struct group *lastg;
  double H;
  int iter, nIter, step, nStep, decorStep;
  double **predA, Z=0.0;
  int i, j;
  int iterFactor = 250;
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
    DECORRELATION STEP
  */
  partRef = CopyPartition(part);
  decorStep = 0;
  fprintf(stderr, "%d 1.0\n", decorStep);
  do {
    for (step=0; step<nnod; step++) {
      MissingLinksMCStep(&H, nlist, glist, part,
			 nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
			 gen);
      decorStep += 1;
    }
    mutualInfo = MutualInformation(partRef, part);
    fprintf(stderr, "%d %g\n", decorStep, mutualInfo);
  } while (mutualInfo > 0.5 && decorStep < 1e5);

/*   /\* */
/*     THERMALIZATION */
/*   *\/ */
/*   H = PartitionH(part); */
/*   nStep = nnod * (int)sqrt(nnod); */
/*   for (iter=0; iter<10; iter++) { */
/*     for (step=0; step<nStep; step++) { */
/*       MissingLinksMCStep(&H, nlist, glist, part, */
/* 			 nnod, G2G, n2gList, LogChooseList, LogChooseListSize, */
/* 			 gen); */
/*     } */
/*     fprintf(stderr, "%d %lf\n", iter, H); */
/*   } */


/*   /\* */
/*     SAMPLING */
/*   *\/ */
/*   nIter = iterFactor * nnod; */
/*   nStep = nnod; */
/*   for (iter=0; iter<nIter; iter++) { */
/*     for (step=0; step<nStep; step++) { */
/*       MissingLinksMCStep(&H, nlist, glist, part, */
/* 			 nnod, G2G, n2gList, LogChooseList, LogChooseListSize, */
/* 			 gen); */
/*     } */

/*     if (iter % 100 == 0) */
/*       fprintf(stderr, "%d %lf\n", iter, H); */
/* /\*       fprintf(stderr, "%d %lf %lf\n", iter, H, PartitionH(part)); *\/ */

/*     /\* Update partition function *\/ */
/*     weight = exp(-H); */
/*     Z += weight; */

/*     /\* Update the predicted adjacency matrix by going through all */
/*        group pairs *\/ */
/*     for (i=0; i<nnod; i++) { */
/*       if (glist[i]->size > 0) { */
	
/* 	/\* update the within-group pairs *\/ */
/* 	r = glist[i]->size * (glist[i]->size - 1) / 2; */
/* 	l = glist[i]->inlinks; */
/* 	contrib = weight * (float)(l + 1) / (float)(r + 2); */
/* 	p1 = glist[i]->nodeList; */
/* 	while ((p1 = p1->next) != NULL) { */
/* 	  p2 = p1; */
/* 	  while ((p2 = p2->next) != NULL) { */
/* 	    predA[p1->node][p2->node] += contrib; */
/* 	    predA[p2->node][p1->node] += contrib; */
/* 	  } */
/* 	} */
      
/* 	/\* update the between-group pairs *\/ */
/* 	for (j=i+1; j<nnod; j++) { */
/* 	  if (glist[j]->size > 0) { */
/* 	    l = G2G[i][j]; */
/* 	    r = glist[i]->size * glist[j]->size; */
/* 	    contrib = weight * (float)(l + 1) / (float)(r + 2); */
/* 	    p1 = glist[i]->nodeList; */
/* 	    while ((p1 = p1->next) != NULL) { */
/* 	      p2 = glist[j]->nodeList; */
/* 	      while ((p2 = p2->next) != NULL) { */
/* 		predA[p1->node][p2->node] += contrib; */
/* 		predA[p2->node][p1->node] += contrib; */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */
/*     } /\* Done updating adjacency matrix *\/ */


/* /\*     for (i=0; i<nnod; i++) { *\/ */
/* /\*       for (j=i+1; j<nnod; j++) { *\/ */
/* /\* 	if (nlist[i]->inGroup == nlist[j]->inGroup) { *\/ */
/* /\* 	  l = glist[nlist[i]->inGroup]->inlinks; *\/ */
/* /\* 	  r = glist[nlist[i]->inGroup]->size * *\/ */
/* /\* 	    (glist[nlist[i]->inGroup]->size - 1) / 2; *\/ */
/* /\* 	} *\/ */
/* /\* 	else { *\/ */
/* /\* 	  l = G2G[nlist[i]->inGroup][nlist[j]->inGroup]; *\/ */
/* /\* 	  r = glist[nlist[i]->inGroup]->size * glist[nlist[j]->inGroup]->size; *\/ */
/* /\* 	} *\/ */
/* /\* 	predA[i][j] += weight * (float)(l + 1) / (float)(r + 2); *\/ */
/* /\* 	predA[j][i] += weight * (float)(l + 1) / (float)(r + 2); *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */

/*   }  /\* End of iter loop *\/ */

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
