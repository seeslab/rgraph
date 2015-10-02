#include <gsl/gsl_rng.h>
#include "partition.h"

unsigned int
GeneralSA(Partition **ppart, AdjaArray *adj,
		  double fac,
		  double Ti, double Tf, double Ts,
		  double cluster_prob,
		  unsigned int nochange_limit,
		  gsl_rng *gen);

unsigned int
SplitModuleSA(unsigned int target, unsigned int empty,
			  double Ti, double Tf, double Ts,
			  unsigned int nochange_limit,
			  Partition *part, AdjaArray *adj,
			  gsl_rng *gen);
