/*
  multiblock.h
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
double PartitionHMB(struct group *part, double linC, double *HarmonicList);
void LSMCStepMB(int factor,
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
int GetDecorrelationStepMB(double *H,
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
void ThermalizeLSMCMB(int decorStep,
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
double **LinkScoreMB(struct node_gra *net,
		      double linC,
		      int nIter,
		      gsl_rng *gen,
		      char verbose_sw);
double ORPartitionHMB(struct group *part, double linC, double *HarmonicList);
void ORLSMCStepMB(int factor,
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
int ORGetDecorrelationStepMB(double *H,
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
void ORThermalizeLSMCMB(int decorStep,
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
double **ORLinkScoreMB(struct node_gra *net,
                      double linC,
                      int nIter,
                      gsl_rng *gen,
                      char verbose_sw);

/*
  ---------------------------------------------------------------------
  Link reliability Gibbs sampling
  ---------------------------------------------------------------------
*/
void GibbsLinkScoreStepMB(double *H,
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
void GibbsThermalizeLinkScoreMB(double *H,
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
double **GibbsLinkScoreMB(struct node_gra *net,
			double linC,
			int nIter,
			gsl_rng *gen,
			char verbose_sw);
void ORGibbsLinkScoreStepMB(double *H,
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
void ORGibbsThermalizeLinkScoreMB(double *H,
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
double **ORGibbsLinkScoreMB(struct node_gra *net,
			double linC,
			int nIter,
			gsl_rng *gen,
			char verbose_sw);


#endif /* !RGRAPH_MULTIBLOCK_H */
