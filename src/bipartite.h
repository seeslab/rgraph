/*
  bipartite.h
  $LastChangedDate: 2007-06-26 16:30:04 -0500 (Tue, 26 Jun 2007) $
  $Revision: 52 $
*/

#ifndef RGRAPH_BIPARTITE_H
#define RGRAPH_BIPARTITE_H 1

#include "prng.h"
#include "graph.h"

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
struct binet *CreateBinet();
struct binet *FBuildNetworkBipart(FILE *inFile,
				  int weight_sw,
				  int add_weight_sw);
void RemoveBinet(struct binet *net);

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
struct binet *InvertBinet(struct binet *net);
struct binet *CopyBinet(struct binet *binet);
struct node_gra *ProjectBinet(struct binet *binet);



#endif /* !RGRAPH_BIPARTITE_H */
