/*
  models.h
  $LastChangedDate$
  $Revision$
*/

#ifndef RGRAPH_MODELS_H
#define RGRAPH_MODELS_H 1

#include "prng.h"

struct node_gra *EmptyGraph(int S);
struct node_gra *ERGraph(int S, double p, struct prng *gen);
struct node_gra *PAGraph(int S, int m, struct prng *gen);
struct node_gra *UndirectedBlockGraph(int ngroup,
				      int *gsize,
				      double **q,
				      struct prng *gen);
struct node_gra *GirvanNewmanGraph(int ngroup,
				   int gsize,
				   double kin,
				   double kout,
				   struct prng *gen);

#endif /* !RGRAPH_MODELS_H */
