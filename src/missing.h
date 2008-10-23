/*
  missing.h
  $LastChangedDate: 2008-10-13 19:13:23 -0500 (Mon, 13 Oct 2008) $
  $Revision: 130 $
*/

#ifndef RGRAPH_MISSING_H
#define RGRAPH_MISSING_H 1

#include "prng.h"

/*
  ---------------------------------------------------------------------
  Missing links
  ---------------------------------------------------------------------
*/
double PartitionH(struct group *part);
void
MissingLinksMCStep(double *H,
		   struct node_gra **nlist,
		   struct group **glist,
		   struct group *part,
		   int nnod,
		   int **G2G,
		   int *n2gList,
		   double **LogChooseList,
		   int LogChooseListSize,
		   struct prng *gen);
double **MissingLinks(struct node_gra *net, struct prng *gen);







#endif /* !RGRAPH_MISSING_H */
