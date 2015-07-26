#include <gsl/gsl_rng.h>
#include "graph.h"
#include "modules.h"
#include "bipartite.h"
#include "partition.h"

struct group *
NetworkClustering(struct node_gra *net,
				  double fac,
				  double Ti, double Tf, double Ts,
				  double cluster_prob,
				  unsigned int nochange_limit, unsigned int Ngroups,
				  gsl_rng *gen);

struct group *
BipartiteNetworkClustering(struct binet *binet,
						   double fac,
						   double Ti, double Tf, double Ts,
						   double cluster_prob,
						   unsigned int nochange_limit, unsigned int Ngroups,
						   unsigned int weighted,
						   gsl_rng *gen);

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




