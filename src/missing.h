/*
  missing.h
  $LastChangedDate: 2008-10-13 19:13:23 -0500 (Mon, 13 Oct 2008) $
  $Revision: 130 $
*/

#ifndef RGRAPH_MISSING_H
#define RGRAPH_MISSING_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include <gsl/gsl_rng.h>
#include "graph.h"
#include "modules.h"

/*
  ---------------------------------------------------------------------
  Missing links
  ---------------------------------------------------------------------
*/
int ExponentialRootF(const gsl_vector *params, void *points, gsl_vector *f);
double CalculateDecay(int nnod, double x1, double y1, double x2, double y2);
double PartitionH(struct group *part, double linC);
void LinkScoreMCStep(double *H,
		     double linC,
		     struct node_gra **nlist,
		     struct group **glist,
		     struct group *part,
		     int nnod,
		     int **G2G,
		     int *n2gList,
		     double **LogChooseList,
		     int LogChooseListSize,
		     gsl_rng *gen);
int GetDecorrelationStep(double *H,
			 double linC,
			 struct node_gra **nlist,
			 struct group **glist,
			 struct group *part,
			 int nnod,
			 int **G2G,
			 int *n2gList,
			 double **LogChooseList,
			 int LogChooseListSize,
			 gsl_rng *gen,
			 char verbose_sw);
void ThermalizeLinkScoreMC(int decorStep,
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
			   gsl_rng *gen,
			   char verbose_sw);
double **LinkScore(struct node_gra *net,
		   double linC,
		   int nIter,
		   gsl_rng *gen,
		   char verbose_sw);
double SBMError(struct node_gra *net, gsl_rng *gen);
double SBMStructureScore(struct node_gra *net, int nrep, gsl_rng *gen);
struct node_gra *NetFromSBMScores(struct node_gra *net, gsl_rng *gen);
double NetworkScore(struct node_gra *netTar,
		    struct node_gra *netObs,
		    double linC,
		    int nIter,
		    gsl_rng *gen,
		    char verbose_sw);
struct node_gra *NetReconstruct(struct node_gra *netObs,
				gsl_rng *gen);







#endif /* !RGRAPH_MISSING_H */
