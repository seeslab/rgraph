/*
  bipartite.h
  $LastChangedDate$
  $Revision$
*/

#ifndef RGRAPH_BIPARTITE_H
#define RGRAPH_BIPARTITE_H 1

#include <gsl/gsl_rng.h>
#include "graph.h"
#include "modules.h"

/*
  ---------------------------------------------------------------------
  Definition of the binet structure A bipartite network is actually a
  couple of networks whose links are only to the other network
  ---------------------------------------------------------------------
*/
struct binet {
  struct node_gra *net1;   // subnetwork 1
  struct node_gra *net2;   // subnetwork 2
};

/*
  ---------------------------------------------------------------------
  Network creation and removal
  ---------------------------------------------------------------------
*/
struct binet *CreateBipart();
struct binet *BuildModularBipartiteNetwork(int *mod_size,
					   int nodpermod,
					   int nmod,
					   double *col_prob,
					   int S2,
					   int mmin, int mmax,
					   double geom_p,
					   double p,
					   gsl_rng *gen);
struct binet *FBuildNetworkBipart(FILE *inFile,
				  int weight_sw,
				  int add_weight_sw);
void RemoveBipart(struct binet *net);

/*
  ---------------------------------------------------------------------
  Node and multinode operations
  ---------------------------------------------------------------------
*/
int NCommonLinksBipart(struct node_gra *n1, struct node_gra *n2);
double SumProductsOfCommonWeightsBipart(struct node_gra *n1, struct node_gra *n2);
void RemoveNodeBipart(struct binet *binet, char *label, int set);

/*
  ---------------------------------------------------------------------
  Network operations
  ---------------------------------------------------------------------
*/
struct binet *CopyBipart(struct binet *binet);
struct binet *InvertBipart(struct binet *net);
int NLinksBipart(struct binet *binet);
struct node_gra *ProjectBipart(struct binet *binet);
struct node_gra *ProjectBipartWeighted(struct binet *binet);
struct binet *RandomizeBipart(struct binet *binet,
			     double times, gsl_rng *gen);

/*
  ---------------------------------------------------------------------
  Network output
  ---------------------------------------------------------------------
*/
void FPrintPajekFileBipart(char *fname,
			   struct binet *binet,
			   int coor_sw,
			   int weight_sw);

void FPrintBipart (FILE *outf, struct binet *binet, int weight_sw);
void FPrintTabNodesBipart(FILE *outf, struct binet *network,  struct group *modules, int degree_based, int weighted);


/*
  ---------------------------------------------------------------------
  Network modularity
  ---------------------------------------------------------------------
*/
double ModularityBipart(struct binet *binet, struct group *part);
double ModularityBipartWeighted(struct binet *binet, struct group *part);
double ModularityBipartWeightedFast(struct binet *binet, struct group *part, double **swwmat);
double ParticipationCoefficientBipart(struct node_gra *node);
void StatisticsParticipationCoefficientBipart(struct node_gra *net,
					      double *theMean,
					      double *theStddev,
					      double *theMin,
					      double *theMax);

#endif /* !RGRAPH_BIPARTITE_H */
