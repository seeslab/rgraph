/*
*/

#ifndef RGRAPH_SPARSE_MISSING_H
#define RGRAPH_SPARSE_MISSING_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include <gsl/gsl_rng.h>
#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "recommend.h"
#include "missing.h"

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
double SparsePartitionH(struct group *part, double betaA, double betaB);
void SparseGibbsLinkScoreStep(double *H,
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
			      gsl_rng *gen);
void SparseGibbsThermalizeLinkScore(double *H,
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
				    gsl_rng *gen, char verbose_sw);
double **SparseGibbsLinkScore(struct node_gra *net,
			      int nIter,
			      gsl_rng *gen,
			      char verbose_sw);



#endif /* !RGRAPH_SPARSE_MISSING_H */
