/*
  bipartite.h
  $LastChangedDate$
  $Revision$
*/

#ifndef RGRAPH_BIPARTITE_H
#define RGRAPH_BIPARTITE_H 1

#include "prng.h"
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
struct binet * BuildModularBipartiteNetwork(int *mod_size,
					    int nodpermod,
					    int nmod,
					    double *col_prob,
					    int S2,
					    int mmin, int mmax,
					    double geom_p,
					    double p,
					    struct prng *gen);
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


/*
  ---------------------------------------------------------------------
  Network operations
  ---------------------------------------------------------------------
*/
struct binet *CopyBipart(struct binet *binet);
struct binet *InvertBipart(struct binet *net);
int NLinksBipart(struct binet *binet);
struct node_gra *ProjectBipart(struct binet *binet);
struct binet *RandomizeBipart(struct binet *binet,
			     double times, struct prng *gen);

/*
  ---------------------------------------------------------------------
  Network output
  ---------------------------------------------------------------------
*/
void FPrintPajekFileBipart(char *fname,
			   struct binet *binet,
			   int coor_sw,
			   int weight_sw);


/*
  ---------------------------------------------------------------------
  Network modularity
  ---------------------------------------------------------------------
*/
double ModularityBipart(struct binet *binet, struct group *part);
void SAGroupSplitBipart(struct group *target_g, struct group *empty_g,
			double Ti, double Tf, double Ts,
			double cluster_prob, 
			double **cmat, double msfac,
			struct prng *gen);
struct group *SACommunityIdentBipart(struct binet *binet,
				    double Ti, double Tf, double Ts,
				    double fac,
				    int ngroup,
				    char initial_sw,
				    int collective_sw,
				    char output_sw,
				    struct prng *gen);


#endif /* !RGRAPH_BIPARTITE_H */
