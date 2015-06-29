#include <gsl/gsl_rng.h>
#include "graph.h"
#include "modules.h"

// Spliting group using SA
struct group *
SAGroupSplit(struct group *targ,
			   double Ti, double Tf, double Ts,
			   int cluster_sw,
			   int weighted,
			   gsl_rng *gen);


void
SAGroupSplitBipart(struct group *target_g, struct group *empty_g,
					 double Ti, double Tf, double Ts,
					 double cluster_prob,
					 double **cmat, double msfac,
					 gsl_rng *gen);

void
SAGroupSplitBipartWeighted(struct group *target_g, struct group *empty_g,
							 double Ti, double Tf, double Ts,
							 double cluster_prob,
							 double *strength,
							 double **swwmat, double Wafac,
							 gsl_rng *gen);

// Simmulated annealing
struct group *
SACommunityIdent(struct node_gra *net,
				   double Ti, double Tf, double Ts,
				   double fac,
				   int ngroup,
				   char initial_sw,
				   int collective_sw,
				   char output_sw,
				   gsl_rng *gen);
struct group *
SACommunityIdentWeight(struct node_gra *net,
						 double Ti, double Tf, double Ts,
						 double fac,
						 int merge,
						 gsl_rng *gen);

struct group *
SACommunityIdentBipart(struct binet *binet,
						 double Ti, double Tf, double Ts,
						 double fac,
						 int ngroup,
						 char initial_sw,
						 int collective_sw,
						 char output_sw,
						 gsl_rng *gen);

struct group *
SACommunityIdentBipartWeighted(struct binet *binet,
								 double Ti, double Tf, double Ts,
								 double fac,
								 int ngroup,
								 char initial_sw,
								 int collective_sw,
								 char output_sw,
								 gsl_rng *gen);
