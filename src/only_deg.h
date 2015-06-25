/*
  only_deg.h
  $LastChangedDate: 2008-10-13 19:13:23 -0500 (Mon, 13 Oct 2008) $
  $Revision: 130 $
*/

#ifndef RGRAPH_MULTIBLOCK_H
#define RGRAPH_MULTIBLOCK_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include <gsl/gsl_rng.h>
#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "missing.h"

double LogDegeneracy(int ngroup);
/*
  ---------------------------------------------------------------------
  Missing links
  ---------------------------------------------------------------------
*/
double PartitionHMB_OD(struct group *part, double linC, double *HarmonicList);
void LSMCStepMB_OD(int factor,
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
		gsl_rng *gen);
int GetDecorrelationStepMB_OD(double *H,
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
			   char verbose_sw);
void ThermalizeLSMCMB_OD(int decorStep,
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
		      char verbose_sw);
double **LinkScoreMB_OD(struct node_gra *net,
		      double linC,
		      int nIter,
		      gsl_rng *gen,
		      char verbose_sw);


/*
  ---------------------------------------------------------------------
  Link reliability Gibbs sampling
  ---------------------------------------------------------------------
*/
void GibbsLinkScoreStepMB_OD(double *H,
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
			gsl_rng *gen);
void GibbsThermalizeLinkScoreMB_OD(double *H,
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
			      char verbose_sw);
double **GibbsLinkScoreMB_OD(struct node_gra *net,
			double linC,
			int nIter,
			gsl_rng *gen,
			char verbose_sw);

#endif /* !RGRAPH_ONLY_DEG_H */
