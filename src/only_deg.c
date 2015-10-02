/*
  only_deg.c
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
#include "multiblock.h"

#define EPSILON 1.e-6

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
double LOGDEGENERACY_OD[150] = {0.0, 0.0, 1.09861228866811, 2.7080502011, 4.7273878187, 7.0501225203, 9.6241042829, 12.4123914513, 15.3873302965, 18.5278179713, 21.8171996489, 25.2419286042, 28.7907695492, 32.454252335, 36.2242817678, 40.0938550791, 44.0568528741, 48.1078805122, 52.242145427, 56.4553607847, 60.7436687672, 65.1035787022, 69.5319165992, 74.0257835705, 78.5825212602, 83.199682859, 87.8750086171, 92.6064050134, 97.3919269199, 102.229762241, 107.118218613, 112.055711822, 117.040755682, 122.07195314, 127.147988427, 132.267620118, 137.429674949, 142.633042317, 147.876669353, 153.1595565, 158.480753534, 163.839355975, 169.234501835, 174.665368673, 180.131170913, 185.631157403, 191.164609186, 196.730837457, 202.329181691, 207.959007927, 213.619707184, 219.310694002, 225.0314051, 230.781298124, 236.559850499, 242.366558357, 248.200935543, 254.062512691, 259.950836364, 265.865468252, 271.805984418, 277.771974602, 283.763041565, 289.77880047, 295.818878313, 301.882913378, 307.970554727, 314.081461726, 320.215303595, 326.371758979, 333.838111, 340.085763, 346.354094, 352.642816, 358.951649, 365.280322, 371.628571, 377.996137, 384.382770, 390.788225, 397.212264, 403.654655, 410.115171, 416.593591, 423.089700, 429.603286, 436.134145, 442.682075, 449.246880, 455.828368, 462.426352, 469.040649, 475.671079, 482.317466, 488.979641, 495.657433, 502.350680, 509.059219, 515.782894, 522.521550, 529.275035, 536.043202, 542.825904, 549.623000, 556.434349, 563.259815, 570.099263, 576.952561, 583.819580, 590.700193, 597.594275, 604.501703, 611.422358, 618.356120, 625.302875, 632.262509, 639.234908, 646.219964, 653.217567, 660.227613, 667.249995, 674.284612, 681.331362, 688.390146, 695.460865, 702.543425, 709.637729, 716.743685, 723.861201, 730.990186, 738.130552, 745.282212, 752.445078, 759.619067, 766.804094, 774.000077, 781.206935, 788.424588, 795.652957, 802.891964, 810.141534, 817.401590, 824.672058, 831.952865, 839.243938, 846.545207, 853.856601, 861.178050, 868.509486, 875.850842};
double
LogDegeneracy_OD(int ngroup)
{
  if (ngroup < 150) {
    return LOGDEGENERACY_OD[ngroup];
  }
  else {
    return 1.4681 * ngroup * (log(ngroup) - 1.);
  }
}


/*
  ---------------------------------------------------------------------
  Partition H
  ---------------------------------------------------------------------
*/
double
PartitionHMB_OD(struct group *part, double linC, double *HarmonicList)
{
  struct group *g1=part, *g2;
  int r, l;
  int ng=0;      /* Number of non-empty groups */
  int nnod=0;  /* Number of nodes */
  double H=0.0;

  H += linC * NNonEmptyGroups(part);

  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
      ng++;
      nnod += g1->size;
      r = g1->size * (g1->size - 1) / 2;
      l = g1->inlinks;
      H += log(r + 1) + LogChoose(r, l);
      /*H -= log(HarmonicList[r + 1] - HarmonicList[l]);*/
      g2 = g1;
      while ((g2 = g2->next) != NULL) {
	if (g2->size > 0) {
	  r = g1->size * g2->size;
	  l = NG2GLinks(g1, g2);
	  H += log(r + 1) + LogChoose(r, l);
	  /*H -= log(HarmonicList[r + 1] - HarmonicList[l]);*/
	}
      }
    }
  }

  H -= gsl_sf_lnfact(nnod - ng);
  H -= LogDegeneracy_OD(ng);

  return H;
}

/*
  ---------------------------------------------------------------------
  Do a Monte Carlo step for the prediction of missing links. Each
  "step" involves nnod independent attempts to move a single node.
  ---------------------------------------------------------------------
*/
void
LSMCStepMB_OD(int factor,
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
	   double *HarmonicList,
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
  int effectng=0;

  for (move=0; move<nnod*factor; move++) {

    /* Choose node and destination group */
    dice = floor(gsl_rng_uniform(gen) * (double)nnod);
    node = nlist[dice];
    oldgnum = node->inGroup;
    do {
      newgnum = floor(gsl_rng_uniform(gen) * (double)nnod);
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
    effectng = 0;
    while ((g = g->next) != NULL) {
      if (g->size > 0) {  /* group is not empty */
	effectng++;
	n2gList[g->label] = NLinksToGroup(node, g);
	if (g->label == oldg->label) {
	  /* old conf, oldg-oldg */
	  r = noldg * (noldg - 1) / 2;
	  l = oldg2oldg;
	  dH -= log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /*dH -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* new conf, oldg-oldg */
	  r = (noldg - 1) * (noldg - 2) / 2;
	  l = oldg2oldg - n2oldg;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /*dH += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* old conf, newg-oldg */
	  r = nnewg * noldg;
	  l = newg2oldg;
	  dH -= log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /*dH -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* new conf, newg-oldg */
	  r = (nnewg + 1) * (noldg - 1);
	  l = newg2oldg + n2oldg - n2newg;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /*dH += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	}
	else if (g->label == newg->label) {
	  /* old conf, newg-newg */
	  r = nnewg * (nnewg - 1) / 2;
	  l = newg2newg;
	  dH -= log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /*dH -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* new conf, newg-newg */
	  r = (nnewg + 1) * nnewg / 2;
	  l = newg2newg + n2newg;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /*dH += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
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
	  /*dH -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* new conf, oldg-g */
	  r = (noldg - 1) * ng;
	  l = oldg2g - n2g;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /*dH += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* old conf, newg-g */
	  r = nnewg * ng;
	  l = newg2g;
	  dH -= log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /*dH -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* new conf, newg-g */
	  r = (nnewg + 1) * ng;
	  l = newg2g + n2g;
	  dH += log(r + 1) + FastLogChoose(r, l,
					   LogChooseList, LogChooseListSize);
	  /*dH += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
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
      dH += LogDegeneracy_OD(effectng) - LogDegeneracy_OD(effectng - 1);
    }
    else if (noldg != 1 && nnewg == 0) {  // effectng would increase by one
      dH -= -FastLogFact(ngroup - effectng, LogFactList, LogFactListSize);
      dH += -FastLogFact(ngroup - (effectng+1), LogFactList, LogFactListSize);
      dH += LogDegeneracy_OD(effectng) - LogDegeneracy_OD(effectng + 1);
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
    }
  } /* nnod moves completed: done! */
}


/*
  ---------------------------------------------------------------------
  
  ---------------------------------------------------------------------
*/
int
GetDecorrelationStepMB_OD(double *H,
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
		       double *HarmonicList,
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
      LSMCStepMB_OD(1, H, linC, nlist, glist, part,
		 nnod, G2G, n2gList,
		 LogChooseList, LogChooseListSize,
		 LogFactList, LogFactListSize, HarmonicList,
		 gen);
      if (step == x1)
		y1 = MutualInformation(partRef, part,0);
    }
    y2 = MutualInformation(partRef, part,0);
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
ThermalizeLSMCMB_OD(int decorStep,
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
		 double *HarmonicList,
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
      LSMCStepMB_OD(decorStep, H, linC, nlist, glist, part,
		 nnod, G2G, n2gList,
		 LogChooseList, LogChooseListSize,
		 LogFactList, LogFactListSize, HarmonicList,
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
LinkScoreMB_OD(struct node_gra *net,
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
  double *HarmonicList=NULL;
  struct node_lis *p1=NULL, *p2=NULL;
  double contrib;
  int dice;
  double Z=0.0;
  int r, l;
  double mutualInfo;
  int LogFactListSize = 10000;
  double *LogFactList=InitializeFastLogFact(LogFactListSize);
  
  /*
    PRELIMINARIES
  */
  /* Initialize the table of harmonic numbers */
  HarmonicList = InitializeHarmonicList(1 + nnod * (nnod - 1) / 2);

  /* Initialize the predicted adjacency matrix */
  predA = allocate_d_mat(nnod, nnod);
  for (i=0; i<nnod; i++) {
    for (j=0; j<nnod; j++) {
      predA[i][j] = 0.0;
    }
  }

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
  H = PartitionHMB_OD(part, linC, HarmonicList);
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "# CALCULATING DECORRELATION TIME\n");
    fprintf(stderr, "# ------------------------------\n");
    break;
  }
  decorStep = GetDecorrelationStepMB_OD(&H, linC, nlist, glist, part,
  				     nnod, G2G, n2gList,
  				     LogChooseList, LogChooseListSize,
  				     LogFactList, LogFactListSize,
  				     HarmonicList,
  				     gen, verbose_sw);

  /* Thermalization */
  switch (verbose_sw) {
  case 'q':
    break;
  default:
    fprintf(stderr, "#\n#\n# THERMALIZING\n");
    fprintf(stderr, "# ------------\n");
    break;
  }
  ThermalizeLSMCMB_OD(decorStep, &H, linC, nlist, glist, part,
  		   nnod, G2G, n2gList,
  		   LogChooseList, LogChooseListSize,
  		   LogFactList, LogFactListSize, HarmonicList,
  		   gen, verbose_sw);
  
  /*
    SAMPLIN' ALONG
  */
  /* Unless we are in debug mode, reset the origin of energies */
  switch (verbose_sw) {
  case 'd':
    break;
  default:
    H = 0.0;
    break;
  }

  /* Do the MC Steps */
  for (iter=0; iter<nIter; iter++) {
    LSMCStepMB_OD(decorStep, &H, linC, nlist, glist, part,
	       nnod, G2G, n2gList,
	       LogChooseList, LogChooseListSize,
	       LogFactList, LogFactListSize, HarmonicList,
	       gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf\n",
	      iter, H, PartitionHMB_OD(part, linC, HarmonicList));
      break;
    }

    /* Update partition function(s) */
    Z += 1.0;
    

    /* Update the predicted adjacency matrix by going through all
       group pairs */
    for (i=0; i<nnod; i++) {
      if (glist[i]->size > 0) {
	
	/* update the within-group pairs */
	r = glist[i]->size * (glist[i]->size - 1) / 2;
	l = glist[i]->inlinks;
	/*contrib =							\
	  (float)(l + 1) * (HarmonicList[r+2] - HarmonicList[l+1])	\
	  / ((float)(r + 2) * (HarmonicList[r+1] - HarmonicList[l]));*/
	contrib = ((float)(l + 1)) / ((float)(r + 2));
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
	    /*contrib =						\
	      (float)(l + 1) * (HarmonicList[r+2] - HarmonicList[l+1])	\
	      / ((float)(r + 2) * (HarmonicList[r+1] - HarmonicList[l]));*/
	    contrib = ((float)(l + 1)) / ((float)(r + 2));

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
  RemovePartition(part);
  free(glist);
  free(nlist);
  free_i_mat(G2G, nnod);
  free_i_vec(n2gList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  FreeFastLogFact(LogFactList);
  FreeHarmonicList(HarmonicList);

  return predA;
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
GibbsLinkScoreStepMB_OD(double *H,
		   double linC,
		   struct node_gra **nlist,
		   struct group **glist,
		   struct group *part,
		   int nnod,
		   int **G2G,
		   int *n2gList,
		   double **LogChooseList,
		   int LogChooseListSize,
		   double *LogFactList, int LogFactListSize, double *HarmonicList,	
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
  double Href;
  
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
	      /*dH[newgnum] -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	      /* new conf, oldg-olg */
	      r = (noldg - 1) * (noldg - 2) / 2;
	      l = oldg2oldg - n2oldg;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /*dH[newgnum] += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	      /* old conf, newg-oldg */
	      r = nnewg * noldg;
	      l = newg2oldg;
	      dH[newgnum] -= log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /*dH[newgnum] -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	      /* new conf, newg-oldg */
	      r = (nnewg + 1) * (noldg - 1);
	      l = newg2oldg + n2oldg - n2newg;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /*dH[newgnum] += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	    }
	    else if (g->label == newg->label) {
	      /* old conf, newg-newg */
	      r = nnewg * (nnewg - 1) / 2;
	      l = newg2newg;
	      dH[newgnum] -= log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /*dH[newgnum] -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	      /* new conf, newg-olg */
	      r = (nnewg + 1) * nnewg / 2;
	      l = newg2newg + n2newg;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /*dH[newgnum] += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
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
	      /*dH[newgnum] -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	      /* new conf, oldg-g */
	      r = (noldg - 1) * ng;
	      l = oldg2g - n2g;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /*dH[newgnum] += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	      /* old conf, newg-g */
	      r = nnewg * ng;
	      l = newg2g;
	      dH[newgnum] -= log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /*dH[newgnum] -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	      /* new conf, newg-g */
	      r = (nnewg + 1) * ng;
	      l = newg2g + n2g;
	      dH[newgnum] += log(r + 1) + FastLogChoose(r, l,
							LogChooseList,
							LogChooseListSize);
	      /*dH[newgnum] += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	    }
	  }
	  else { /* group is empty */
	    n2gList[g->label] = 0.0;
	  }
	}


	/* Labeled-groups sampling correction */
	if (noldg == 1 && nnewg != 0) { // effectng would decrease by one
	  dH[newgnum] -= -FastLogFact(ngroup - effectng, LogFactList, LogFactListSize);
	  dH[newgnum] += -FastLogFact(ngroup - (effectng-1), LogFactList, LogFactListSize);
	  dH[newgnum] += LogDegeneracy_OD(effectng) - LogDegeneracy_OD(effectng - 1);
	}
	else if (noldg != 1 && nnewg == 0) {  // effectng would increase by one
	  dH[newgnum] -= -FastLogFact(ngroup - effectng, LogFactList, LogFactListSize);
	  dH[newgnum] += -FastLogFact(ngroup - (effectng+1), LogFactList, LogFactListSize);
	  dH[newgnum] += LogDegeneracy_OD(effectng) - LogDegeneracy_OD(effectng + 1);
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
	  /*dHempty -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* new conf, oldg-olg */
	  r = (noldg - 1) * (noldg - 2) / 2;
	  l = oldg2oldg - n2oldg;
	  dHempty += log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /*dHempty += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* old conf, newg-oldg */
	  r = nnewg * noldg;
	  l = newg2oldg;
	  dHempty -= log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /*dHempty -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* new conf, newg-oldg */
	  r = (nnewg + 1) * (noldg - 1);
	  l = newg2oldg + n2oldg - n2newg;
	  dHempty += log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /*dHempty += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
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
	  /*dHempty -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* new conf, oldg-g */
	  r = (noldg - 1) * ng;
	  l = oldg2g - n2g;
	  dHempty += log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /*dHempty += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* old conf, newg-g */
	  r = nnewg * ng;
	  l = newg2g;
	  dHempty -= log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /*dHempty -= -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	  /* new conf, newg-g */
	  r = (nnewg + 1) * ng;
	  l = newg2g + n2g;
	  dHempty += log(r + 1) + FastLogChoose(r, l,
						LogChooseList,
						LogChooseListSize);
	  /*dHempty += -log(HarmonicList[r + 1] - HarmonicList[l]);*/
	}
      }
    }

    /* Labeled-groups sampling correction */
    if (noldg == 1 && nnewg != 0) { // effectng would decrease by one
      dHempty -= -FastLogFact(ngroup - effectng, LogFactList, LogFactListSize);
      dHempty += -FastLogFact(ngroup - (effectng-1), LogFactList, LogFactListSize);
      dHempty += LogDegeneracy_OD(effectng) - LogDegeneracy_OD(effectng - 1);
    }
    else if (noldg != 1 && nnewg == 0) {  // effectng would increase by one
      dHempty -= -FastLogFact(ngroup - effectng, LogFactList, LogFactListSize);
      dHempty += -FastLogFact(ngroup - (effectng+1), LogFactList, LogFactListSize);
      dHempty += LogDegeneracy_OD(effectng) - LogDegeneracy_OD(effectng + 1);
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
GibbsThermalizeLSMCMB_OD(double *H,
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
		         double *HarmonicList,	
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
      GibbsLinkScoreStepMB_OD(H, linC, nlist, glist, part,
			 nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
			 LogFactList, LogFactListSize, HarmonicList,
			 gen);
      switch (verbose_sw) {
      case 'q':
	/*PartitionHMB(part, linC, HarmonicList);*/ /* Number of groups */
	break;
      case 'd':
	fprintf(stderr, "%lf %lf %d\n", *H, PartitionHMB(part, linC, HarmonicList),
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
GibbsLinkScoreMB_OD(struct node_gra *net,
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
  double Z=0.0;
  int r, l;
  double mutualInfo;
  int LogFactListSize = 10000;
  double *LogFactList = InitializeFastLogFact(LogFactListSize);
  double *HarmonicList = InitializeHarmonicList(1 + nnod * (nnod - 1) / 2);

  /*
    PRELIMINARIES
  */

  /* Initialize the predicted adjacency matrix */
  predA = allocate_d_mat(nnod, nnod);
  for (i=0; i<nnod; i++) {
    for (j=0; j<nnod; j++) {
      predA[i][j] = 0.0;
    }
  }

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
  H = PartitionHMB(part, linC, HarmonicList);

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
  GibbsThermalizeLSMCMB_OD(&H, linC, nlist, glist, part,
			   nnod, G2G, n2gList,
			   LogChooseList, LogChooseListSize,
			   LogFactList, LogFactListSize, HarmonicList,
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
    GibbsLinkScoreStepMB_OD(&H, linC, nlist, glist, part,
		       nnod, G2G, n2gList, LogChooseList, LogChooseListSize,
		       LogFactList, LogFactListSize, HarmonicList,
		       gen);
    switch (verbose_sw) {
    case 'q':
      break;
    case 'v':
      fprintf(stderr, "%d %lf\n", iter, H);
      break;
    case 'd':
      fprintf(stderr, "%d %lf %lf %d\n", iter, H, PartitionHMB(part, linC, HarmonicList),
	      NNonEmptyGroups(part));
      /* FPrintPartition(stderr, part, 0); */
      break;
    }


    /* Update partition function(s) */
    Z += 1.0;


    /* Update the predicted adjacency matrix by going through all
       group pairs */
    for (i=0; i<nnod; i++) {
      if (glist[i]->size > 0) {
	
	/* update the within-group pairs */
	r = glist[i]->size * (glist[i]->size - 1) / 2;
	l = glist[i]->inlinks;
	/*contrib =							\
	  (float)(l + 1) * (HarmonicList[r+2] - HarmonicList[l+1])	\
	  / ((float)(r + 2) * (HarmonicList[r+1] - HarmonicList[l]));*/

	contrib = ((float)(l + 1)) / ((float)(r + 2));
	
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
	    /*contrib =						\
	      (float)(l + 1) * (HarmonicList[r+2] - HarmonicList[l+1])	\
	      / ((float)(r + 2) * (HarmonicList[r+1] - HarmonicList[l]));*/
	    contrib = ((float)(l + 1)) / ((float)(r + 2));

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
  RemovePartition(part);
  free(glist);
  free(nlist);
  free_i_mat(G2G, nnod);
  free_i_vec(n2gList);
  FreeFastLogChoose(LogChooseList, LogChooseListSize);
  FreeFastLogFact(LogFactList);
  FreeHarmonicList(HarmonicList);

  return predA;
}
