/*
  conflict.h
  $LastChangedDate: 2008-10-13 19:13:23 -0500 (Mon, 13 Oct 2008) $
  $Revision: 130 $
*/

#ifndef RGRAPH_CONFLICT_H
#define RGRAPH_CONFLICT_H 1

#include <gsl/gsl_rng.h>
#include "graph.h"
#include "modules.h"
#include "bipartite.h"
#include "recommend.h"

struct query **AllLinkScore2State(struct binet *ratings,
				  int nIter,
				  gsl_rng *gen,
				  char verbose_sw,
				  int decorStep);

#endif /* !RGRAPH_CONFLICT_H */
