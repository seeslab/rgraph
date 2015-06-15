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
#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "recommend.h"

#define EPSILON 1.e-6

/*
  ---------------------------------------------------------------------
  Missing links
  ---------------------------------------------------------------------
*/
int ExponentialRootF(const gsl_vector *params, void *points, gsl_vector *f);
double CalculateDecay(int nnod, double x1, double y1, double x2, double y2);
double PartitionH(struct group *part, double linC);
void LinkScoreMCStep(int factor,
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
			 double *LogFactList, int LogFactListSize,
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
			   double *LogFactList, int LogFactListSize,
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
struct group **PartitionSampling(struct node_gra *net,
				 double linC,
				 int nIter,
				 gsl_rng *gen,
				 char verbose_sw,
				 int burnin,
				 double thinning);

/*
  ---------------------------------------------------------------------
  Link reliability K-state
  ---------------------------------------------------------------------
*/
double LSHKState(int K, struct group *part);
void LSMCStepKState(int K,
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
		    gsl_rng *gen);
int LSGetDecorrelationStepKState(int K,
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
				 char verbose_sw);
void LSThermalizeMCKState(int K,
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
			  double *LogFactList,int LogFactListSize,
			  gsl_rng *gen,
			  char verbose_sw);
double **LSMultiLinkScoreKState(int K,
				struct node_gra *ratings,
				int nIter,
				gsl_rng *gen,
				char verbose_sw,
				int decorStep,
				struct query ***queries,
				int *nquery);



/*
  ---------------------------------------------------------------------
  Link reliability Gibbs sampling
  ---------------------------------------------------------------------
*/
void GibbsLinkScoreStep(double *H,
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
			gsl_rng *gen);
void GibbsThermalizeLinkScore(double *H,
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
			      char verbose_sw);
double **GibbsLinkScore(struct node_gra *net,
			double linC,
			int nIter,
			gsl_rng *gen,
			char verbose_sw);



#endif /* !RGRAPH_MISSING_H */
