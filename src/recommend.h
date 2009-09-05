/*
  recommend.h
  $LastChangedDate: 2008-10-13 19:13:23 -0500 (Mon, 13 Oct 2008) $
  $Revision: 130 $
*/

#ifndef RGRAPH_RECOMMEND_H
#define RGRAPH_RECOMMEND_H 1

#include "prng.h"
#include "graph.h"
#include "modules.h"
#include "bipartite.h"

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
  
  -----------------------------------------------------------------------------
*/
  double H2State(struct group *part1, struct group *part2,
	       struct query **query_list, int nquery);
void MCStep2State(int factor,
		  double *H,
		  struct query **query_list, int nquery, 
		  struct node_gra **nlist1, struct node_gra **nlist2,
		  struct group **glist1, struct group **glist2,
		  struct group *part1, struct group *part2,
		  int nnod1, int nnod2,
		  int **G1G2, int **G2G1,
		  int *n2gList,
		  double **LogChooseList,
		  int LogChooseListSize,
		  struct prng *gen);
double LinkScore2State(struct binet *binet,
		       struct query *the_query,
		       int nIter,
		       struct prng *gen,
		       char verbose_sw);


#endif /* !RGRAPH_RECOMMEND_H */
