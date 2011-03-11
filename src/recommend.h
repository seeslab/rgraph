/*
  recommend.h
  $LastChangedDate: 2008-10-13 19:13:23 -0500 (Mon, 13 Oct 2008) $
  $Revision: 130 $
*/

#ifndef RGRAPH_RECOMMEND_H
#define RGRAPH_RECOMMEND_H 1

#include <gsl/gsl_rng.h>
#include "graph.h"
#include "modules.h"
#include "bipartite.h"

#define MAX_QUERIES 10000;

/*
  -----------------------------------------------------------------------------
  A structure containing a "query link", i.e. an unobserved link for
  which we want to make a recommendation
  -----------------------------------------------------------------------------
*/
struct query
{
  struct node_gra *n1;
  struct node_gra *n2;
};


/*
  -----------------------------------------------------------------------------
  Auxiliary functions
  -----------------------------------------------------------------------------
*/
struct query *CreateQuery(struct node_gra *node1, struct node_gra *node2);
void FreeQuery(struct query *q);
int CountUnobserved(struct binet *ratings);
struct query **BuildUnobservedSet(struct binet *ratings);
void RemoveRatings(struct binet *ratings, int r);

/*
  -----------------------------------------------------------------------------
  I/O functions
  -----------------------------------------------------------------------------
*/
struct binet *ReadRecommenderObservations(FILE *inFile);
struct query **ReadQueries(FILE *inFile, int nQueries, struct binet *binet);

/*
  -----------------------------------------------------------------------------
  Recommender functions
  -----------------------------------------------------------------------------
*/
double H2State(struct group *part1, struct group *part2,
	       struct query **ignore_list, int nignore);
void MCStep2State(int factor,
		  double *H,
		  struct query **ignore_list, int nignore, 
		  struct node_gra **nlist1, struct node_gra **nlist2,
		  struct group **glist1, struct group **glist2,
		  struct group *part1, struct group *part2,
		  int nnod1, int nnod2,
		  int **G1G2, int **G2G1,
		  int *n2gList,
		  double **LogChooseList,
		  int LogChooseListSize,
		  gsl_rng *gen);
int GetDecorrelationStep2State(double *H,
			       struct query **ignore_list, int nignore, 
			       struct node_gra **nlist1,
			       struct node_gra **nlist2,
			       struct group **glist1, struct group **glist2,
			       struct group *part1, struct group *part2,
			       int nnod1, int nnod2,
			       int **G1G2, int **G2G1,
			       int *n2gList,
			       double **LogChooseList,
			       int LogChooseListSize,
			       gsl_rng *gen,
			       char verbose_sw);
void ThermalizeMC2State(int decorStep,
			double *H,
			struct query **ignore_list, int nignore, 
			struct node_gra **nlist1, struct node_gra **nlist2,
			struct group **glist1, struct group **glist2,
			struct group *part1, struct group *part2,
			int nnod1, int nnod2,
			int **G1G2, int **G2G1,
			int *n2gList,
			double **LogChooseList,
			int LogChooseListSize,
			gsl_rng *gen,
			char verbose_sw);
double LinkScore2State(struct binet *binet,
		       struct query *the_query,
		       int nIter,
		       gsl_rng *gen,
		       char verbose_sw,
		       int decorStep);
int IsQueryInSet(struct query *q, struct query **set, int nset);
double *MultiLinkScore2State(struct binet *ratings,
			     struct query **querySet, int nquery,
			     int nIter,
			     gsl_rng *gen,
			     char verbose_sw,
			     int decorStep);


double H2StateFast(struct group *part1, struct group *part2);
void MCStep2StateFast(int factor,
		      double *H,
		      struct node_gra **nlist1, struct node_gra **nlist2,
		      struct group **glist1, struct group **glist2,
		      struct group *part1, struct group *part2,
		      int nnod1, int nnod2,
		      int **N1G2_0, int **N2G1_0,
		      int **N1G2_1, int **N2G1_1,
		      int **G1G2_0, int **G2G1_0,
		      int **G1G2_1, int **G2G1_1,
		      double *LogList,
		      int LogListSize,
		      double **LogChooseList,
		      int LogChooseListSize,
		      gsl_rng *gen);
int GetDecorrelationStep2StateFast(double *H,
				   struct node_gra **nlist1,
				   struct node_gra **nlist2,
				   struct group **glist1, struct group **glist2,
				   struct group *part1, struct group *part2,
				   int nnod1, int nnod2,
				   int **N1G2_0, int **N2G1_0,
				   int **N1G2_1, int **N2G1_1,
				   int **G1G2_0, int **G2G1_0,
				   int **G1G2_1, int **G2G1_1,
				   double *LogList,
				   int LogListSize,
				   double **LogChooseList,
				   int LogChooseListSize,
				   gsl_rng *gen,
				   char verbose_sw);
void ThermalizeMC2StateFast(int decorStep,
			    double *H,
			    struct node_gra **nlist1, struct node_gra **nlist2,
			    struct group **glist1, struct group **glist2,
			    struct group *part1, struct group *part2,
			    int nnod1, int nnod2,
			    int **N1G2_0, int **N2G1_0,
			    int **N1G2_1, int **N2G1_1,
			    int **G1G2_0, int **G2G1_0,
			    int **G1G2_1, int **G2G1_1,
			    double *LogList,
			    int LogListSize,
			    double **LogChooseList,
			    int LogChooseListSize,
			    gsl_rng *gen,
			char verbose_sw);
double *MultiLinkScore2StateFast(struct binet *ratings,
				 struct query **querySet, int nquery,
				 int nIter,
				 gsl_rng *gen,
				 char verbose_sw,
				 int decorStep);


#endif /* !RGRAPH_RECOMMEND_H */
