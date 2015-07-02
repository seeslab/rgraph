#include <gsl/gsl_rng.h>
#include "graph.h"
#include "modules.h"
#include "bipartite.h"

struct llist
{
  unsigned int value;
  struct llist *next;
};

struct partition{
  struct llist **group;
  unsigned int *module;
  unsigned int *size;
  unsigned int nmod;
  unsigned int nnod;
};

double **CostMatrix(struct node_gra *net);
double **CostMatrixWeighted(struct node_gra *net);
double **CostMatrixBipart(struct binet *binet);
double **CostMatrixBipartWeighted(struct binet *binet);

void iterationNumber(int nnod,
					 double fac,
					 int *individual_movements,
					 int *collective_movements);
struct group *convertPartitionToGroup(struct partition *part, struct node_gra *net);
unsigned int move(unsigned int val,
				  unsigned int old, unsigned int new,
				  struct partition *part);
void merge(unsigned int g1,
		   unsigned int g2,
		   struct partition *part);
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



struct partition * InitializePartition(unsigned int N);
struct partition *CopySimplePartition(struct partition *part);
void FreeSimplePartition(struct partition *part);

struct group *
GeneralSA(double **cmat, unsigned int N,
		  struct node_gra *net,
		  double fac,
		  double Ti, double Tf, double Ts,
		  gsl_rng *gen);
