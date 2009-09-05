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
#include "recommend.h"

#define max(A, B) ((A > B)? A : B)

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
/* struct pair */
/* { */
/*   double x0; */
/*   double y0; */
/*   double x1; */
/*   double y1; */
/* }; */

/* int */
/* ExponentialRootF(const gsl_vector *params, */
/* 		 void *points, */
/* 		 gsl_vector * f) */
/* { */
/*   const double a = gsl_vector_get(params, 0); */
/*   const double b = gsl_vector_get(params, 1); */

/*   double x0 = ((struct pair *)points)->x0; */
/*   double y0 = ((struct pair *)points)->y0; */
/*   double x1 = ((struct pair *)points)->x1; */
/*   double y1 = ((struct pair *)points)->y1; */
  
/*   const double r0 = y0 - a - (1. - a) * exp(-x0 / b); */
/*   const double r1 = y1 - a - (1. - a) * exp(-x1 / b); */
  
/*   gsl_vector_set (f, 0, r0); */
/*   gsl_vector_set (f, 1, r1); */
  
/*   return GSL_SUCCESS; */
/* } */

/* double */
/* CalculateDecay(int nnod, double x1, double y1, double x2, double y2) */
/* { */
/*   const gsl_multiroot_fsolver_type *T; */
/*   gsl_multiroot_fsolver *s; */
/*   int status; */
/*   size_t i, iter = 0; */
/*   const size_t n = 2; */
/*   struct pair p = {x1, y1, x2, y2}; */
/*   gsl_multiroot_function f = {&ExponentialRootF, n, &p}; */
/*   double x_init[2] = {y2, sqrt(nnod)}; */
/*   gsl_vector *x = gsl_vector_alloc(n); */
/*   double result; */
  
/*   for (i=0; i<n; i++) */
/*     gsl_vector_set(x, i, x_init[i]); */
  
/*   T = gsl_multiroot_fsolver_hybrids; */
/*   s = gsl_multiroot_fsolver_alloc (T, n); */
/*   gsl_multiroot_fsolver_set (s, &f, x); */
/*   do */
/*     { */
/*       iter++; */
/*       status = gsl_multiroot_fsolver_iterate (s); */
/*       if (status)   /\* check if solver is stuck *\/ */
/* 	break; */
/*       status =gsl_multiroot_test_residual(s->f, 1e-7); */
/*     } */
/*   while (status == GSL_CONTINUE && iter < 1000); */
  
/* /\*   fprintf(stderr, "# GSL status = %s\n", gsl_strerror(status)); *\/ */
/*   if (strcmp(gsl_strerror(status), "success") != 0) */
/*     result = -1; */
/*   else */
/*     result = gsl_vector_get(s->x, 1); */

/*   gsl_multiroot_fsolver_free(s); */
/*   gsl_vector_free(x); */

/*   return result; */
/* } */

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
  struct group *oldg, *newg, *g;
  int ***G2G, ***G2Ginv;
  int move;
  double set_ratio;
  struct node_gra *node=NULL;
  int dice, r, l;
  int oldgnum, newgnum, q;
  int i, nnod;

  /* Ratio of moves in each of the sets */
  set_ratio = (double)(nnod1*nnod1) / (double)(nnod1*nnod1 + nnod2*nnod2);

  for (move=0; move<(nnod1+nnod2)*factor; move++) {
    /* The move */
    if (prng_get_next(gen) < set_ratio) { /* move in first set */
      G2G = &G1G2;
      G2Ginv = &G2G1;
      nnod = nnod1;
      dice = floor(prng_get_next(gen) * (double)nnod1);
      node = nlist1[dice];
      oldgnum = node->inGroup;
      do {
	newgnum = floor(prng_get_next(gen) * (double)nnod1);
      } while (newgnum == oldgnum);
      oldg = glist1[oldgnum];
      newg = glist1[newgnum];
      g = part1;
    }
    else {                                /* move in second set */
      G2G = &G2G1;
      G2Ginv = &G1G2;
      nnod = nnod2;
      dice = floor(prng_get_next(gen) * (double)nnod2);
      node = nlist2[dice];
      oldgnum = node->inGroup;
      do {
	newgnum = floor(prng_get_next(gen) * (double)nnod2);
      } while (newgnum == oldgnum);
      oldg = glist2[oldgnum];
      newg = glist2[newgnum];
      g = part2;
    }

    /* The change of energy */
    dH = 0.0;
    while ((g=g->next) != NULL) {
      if (g->size > 0) {  /* group is not empty */
      	n2gList[g->label] = NLinksToGroup(node, g);
	/* old configuration, old group */
	r = oldg->size * g->size;
	l = (*G2G)[oldgnum][g->label];
	for (q=0; q<nquery; q++) {
	  if (query_list[q]->n1->inGroup == oldg->label &&
	      query_list[q]->n2->inGroup == g->label) {
	    r--;
	    l -= IsThereLink(query_list[q]->n1, query_list[q]->n2);
	  }
	}
	dH -= log(r + 1) + LogChoose(r, l);
	/* old configuration, new group */
	r = newg->size * g->size;
	l = (*G2G)[newgnum][g->label];
	for (q=0; q<nquery; q++) {
	  if (query_list[q]->n1->inGroup == newg->label &&
	      query_list[q]->n2->inGroup == g->label) {
	    r--;
	    l -= IsThereLink(query_list[q]->n1, query_list[q]->n2);
	  }
	}
	dH -= log(r + 1) + LogChoose(r, l);
	/* new configuration, old group */
	r = (oldg->size - 1) * g->size;
	l = (*G2G)[oldgnum][g->label] - n2gList[g->label];
	for (q=0; q<nquery; q++) {
	  if (query_list[q]->n1->inGroup == oldg->label &&
	      query_list[q]->n2->inGroup == g->label) {
	    r--;
	    l -= IsThereLink(query_list[q]->n1, query_list[q]->n2);
	  }
	}
	dH += log(r + 1) + LogChoose(r, l);
	/* new configuration, new group */
	r = (newg->size + 1) * g->size;
	l = (*G2G)[newgnum][g->label] + n2gList[g->label];
	for (q=0; q<nquery; q++) {
	  if (query_list[q]->n1->inGroup == newg->label &&
	      query_list[q]->n2->inGroup == g->label) {
	    r--;
	    l -= IsThereLink(query_list[q]->n1, query_list[q]->n2);
	  }
	}
	dH += log(r + 1) + LogChoose(r, l);
      }
      else { /* group is empty */
	n2gList[g->label] = 0.0;
      }
    }
    
    /* Metropolis acceptance */
    if ((dH <= 0.0) || (prng_get_next(gen) < exp(-dH))) {
      /* accept move */
      MoveNode(node, oldg, newg);
      *H += dH;
      /* update G2G and G2Ginv*/
      for (i=0; i<nnod; i++) {
	(*G2G)[i][oldgnum] -= n2gList[i];
	(*G2G)[i][newgnum] += n2gList[i];
	(*G2Ginv)[oldgnum][i] -= n2gList[i];
	(*G2Ginv)[newgnum][i] += n2gList[i];
      }
    }
  } /* Moves completed: done! */

  return;
}


/* /\* */
/*   --------------------------------------------------------------------- */
  
/*   --------------------------------------------------------------------- */
/* *\/ */
/* int */
/* GetDecorrelationStep(double *H, */
/* 		     double linC, */
/* 		     struct node_gra **nlist, */
/* 		     struct group **glist, */
/* 		     struct group *part, */
/* 		     int nnod, */
/* 		     int **G2G, */
/* 		     int *n2gList, */
/* 		     double **LogChooseList, */
/* 		     int LogChooseListSize, */
/* 		     struct prng *gen, */
/* 		     char verbose_sw) */
/* { */
/*   struct group *partRef; */
/*   int step, x1, x2; */
/*   double y1, y2; */
/*   double mutualInfo; */
/*   int rep, nrep=10; */
/*   double *decay, meanDecay, sigmaDecay, result; */
/*   int norm=0; */

/*   x2 = nnod / 5; */
/*   x1 = x2 / 4; */

/*   /\* Get the nrep initial estimates *\/ */
/*   decay = allocate_d_vec(nrep); */
/*   for (rep=0; rep<nrep; rep++) { */
/*     switch (verbose_sw) { */
/*     case 'q': */
/*       break; */
/*     default: */
/*       fprintf(stderr, "#\n# Estimating decorrelation time (%d/%d)\n", */
/* 	      rep + 1, nrep); */
/*       break; */
/*     } */
/*     partRef = CopyPartition(part); */
/*     for (step=0; step<=x2; step++) { */
/*       LinkScoreMCStep(1, H, linC, nlist, glist, part, */
/* 			 nnod, G2G, n2gList, LogChooseList, LogChooseListSize, */
/* 			 gen); */
/*       if (step == x1) */
/* 	y1 = MutualInformation(partRef, part); */
/*     } */
/*     y2 = MutualInformation(partRef, part); */
/*     RemovePartition(partRef); */
/*     decay[rep] = 2. * CalculateDecay(nnod, x1, y1, x2, y2); */
/*     switch (verbose_sw) { */
/*     case 'q': */
/*       break; */
/*     default: */
/*       fprintf(stderr, "# Decorrelation time (estimate %d) = %g\n", */
/* 	      rep + 1, decay[rep]); */
/*       break; */
/*     } */
/*     if (decay[rep] < 0) { */
/*       rep--; */
/*       switch (verbose_sw) { */
/*       case 'q': */
/* 	break; */
/*       default: */
/* 	fprintf(stderr, "#\tignoring...\n"); */
/* 	break; */
/*       } */
/*     } */
/*   } */
  
/*   /\* Get rid of bad estimates (Chauvenet criterion)  *\/ */
/*   meanDecay = mean(decay, nrep); */
/*   sigmaDecay = stddev(decay, nrep); */
/*   result = meanDecay * nrep; */
/*   for (rep=0; rep<nrep; rep++) { */
/*     if (fabs(decay[rep] - meanDecay) / sigmaDecay > 2) { */
/*       result -= decay[rep]; */
/*       switch (verbose_sw) { */
/*       case 'q': */
/* 	break; */
/*       default: */
/* 	fprintf(stderr, "# Disregarding estimate %d\n", rep + 1); */
/* 	break; */
/*       } */
/*     } */
/*     else { */
/*       norm++; */
/*     } */
/*   } */
  
/*   /\* Clean up *\/ */
/*   free_d_vec(decay); */

/*   return result / norm; */
/* } */

/* /\* */
/*   --------------------------------------------------------------------- */
  
/*   --------------------------------------------------------------------- */
/* *\/ */
/* void */
/* ThermalizeLinkScoreMC(int decorStep, */
/* 		      double *H, */
/* 		      double linC, */
/* 		      struct node_gra **nlist, */
/* 		      struct group **glist, */
/* 		      struct group *part, */
/* 		      int nnod, */
/* 		      int **G2G, */
/* 		      int *n2gList, */
/* 		      double **LogChooseList, */
/* 		      int LogChooseListSize, */
/* 		      struct prng *gen, */
/* 		      char verbose_sw) */
/* { */
/*   double HMean0=1.e10, HStd0=1.e-10, HMean1, HStd1, *Hvalues; */
/*   int rep, nrep=20; */
/*   int equilibrated=0; */

/*   Hvalues = allocate_d_vec(nrep); */

/*   do { */
    
/*     /\* MC steps *\/ */
/*     for (rep=0; rep<nrep; rep++) { */
/*       LinkScoreMCStep(decorStep, H, linC, nlist, glist, part, */
/* 			 nnod, G2G, n2gList, LogChooseList, LogChooseListSize, */
/* 			 gen); */
/*       switch (verbose_sw) { */
/*       case 'q': */
/* 	break; */
/*       default: */
/* 	fprintf(stderr, "%lf\n", *H); */
/* 	break; */
/*       } */
/*       Hvalues[rep] = *H; */
/*     } */

/*     /\* Check for equilibration *\/ */
/*     HMean1 = mean(Hvalues, nrep); */
/*     HStd1 = stddev(Hvalues, nrep); */
/*     if (HMean0 - HStd0 / sqrt(nrep) < HMean1 + HStd1 / sqrt(nrep)) { */
/*       equilibrated++; */
/*       switch (verbose_sw) { */
/*       case 'q': */
/* 	break; */
/*       default: */
/* 	fprintf(stderr, "#\tequilibrated (%d/5) H=%lf\n", */
/* 		equilibrated, HMean1); */
/* 	break; */
/*       } */
/*     } */
/*     else { */
/*       switch (verbose_sw) { */
/*       case 'q': */
/* 	break; */
/*       default: */
/* 	fprintf(stderr, "#\tnot equilibrated yet H0=%g+-%g H1=%g+-%g\n", */
/* 		HMean0, HStd0 / sqrt(nrep), HMean1, HStd1 / sqrt(nrep)); */
/* 	break; */
/*       } */
/*       HMean0 = HMean1; */
/*       HStd0 = HStd1; */
/*       equilibrated = 0; */
/*     } */

/*   } while (equilibrated < 5); */
  
/*   /\* Clean up *\/ */
/*   free_d_vec(Hvalues); */

/*   return; */
/* } */

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
		char verbose_sw)
{
  int nnod1=CountNodes(binet->net1), nnod2=CountNodes(binet->net2);
  struct node_gra *net1=NULL, *net2=NULL;
  struct group *part1=NULL, *part2=NULL;
  struct node_gra *p1=NULL, *p2=NULL, *node=NULL;
  struct node_gra **nlist1=NULL, **nlist2=NULL;
  struct group **glist1=NULL, **glist2=NULL;
  struct group *lastg=NULL;
  double H;
  int iter, decorStep;
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
  int nquery;

  /*
    PRELIMINARIES
  */
  /* The query */
  query_list[0] = the_query;
  nquery = 1;

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
  n2gList = allocate_i_vec(max(nnod1, nnod2));
  for (i=0; i<nnod1; i++) {
    for (j=0; j<nnod2; j++) {
      G1G2[i][j] = G2G1[j][i] = NG2GLinks(glist1[i], glist2[j]);
    }
  }

  /*
    GET READY FOR THE SAMPLING
  */
/*   /\* Get the decorrelation time *\/ */
/*   H = PartitionH(part); */
/*   switch (verbose_sw) { */
/*   case 'q': */
/*     break; */
/*   default: */
/*     fprintf(stderr, "# CALCULATING DECORRELATION TIME\n"); */
/*     fprintf(stderr, "# ------------------------------\n"); */
/*     break; */
/*   } */
/*   decorStep = GetDecorrelationStep(&H, nlist, glist, part, */
/* 				   nnod, G2G, n2gList, */
/* 				   LogChooseList, LogChooseListSize, */
/* 				   gen, verbose_sw); */

/*   /\* Thermalization *\/ */
/*   switch (verbose_sw) { */
/*   case 'q': */
/*     break; */
/*   default: */
/*     fprintf(stderr, "#\n#\n# THERMALIZING\n"); */
/*     fprintf(stderr, "# ------------\n"); */
/*     break; */
/*   } */
/*   ThermalizeLinkScoreMC(decorStep, &H, nlist, glist, part, */
/* 			nnod, G2G, n2gList, */
/* 			LogChooseList, LogChooseListSize, */
/* 			gen, verbose_sw); */
  
  /*
    SAMPLIN' ALONG
  */
  H = 0; /* Reset the origin of energies to avoid huge exponentials */
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
      break;
    }

    /* Update partition function */
    weight = exp(-H);
    Z += weight;

    /* Update the score */
    l = G1G2[the_query->n1->inGroup][the_query->n2->inGroup] -	\
      IsThereLink(the_query->n1, the_query->n2);
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
  free(query_list);
  
  return score;
}
