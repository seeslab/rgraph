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

#endif /* !RGRAPH_MODELS_H */
