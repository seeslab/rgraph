/*
  missing.h
  $LastChangedDate: 2008-10-13 19:13:23 -0500 (Mon, 13 Oct 2008) $
  $Revision: 130 $
*/

#ifndef RGRAPH_MISSING_H
#define RGRAPH_MISSING_H 1

#include "prng.h"

/*
  ---------------------------------------------------------------------
  Missing links
  ---------------------------------------------------------------------
*/
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
		     struct prng *gen);
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
			 struct prng *gen,
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
			   struct prng *gen,
			   char verbose_sw);
double **LinkScore(struct node_gra *net,
		   double linC,
		   int nIter,
		   struct prng *gen,
		   char verbose_sw);
double SBMError(struct node_gra *net, struct prng *gen);
double SBMStructureScore(struct node_gra *net, int nrep, struct prng *gen);
struct node_gra *NetFromSBMScores(struct node_gra *net, struct prng *gen);
double NetworkScore(struct node_gra *netTar,
		    struct node_gra *netObs,
		    double linC,
		    int nIter,
		    struct prng *gen,
		    char verbose_sw);
struct node_gra *NetReconstruct(struct node_gra *netObs,
				struct prng *gen);







#endif /* !RGRAPH_MISSING_H */
