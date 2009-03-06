/*
  matrix.h
  $LastChangedDate: 2008-10-13 19:13:23 -0500 (Mon, 13 Oct 2008) $
  $Revision: 130 $
*/

#ifndef RGRAPH_MATRIX_H
#define RGRAPH_MATRIX_H 1

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

#include "tools.h"
#include "graph.h"
#include "matrix.h"

int **AdjacencyMat(struct node_gra *net);
int **LaplacianMat(struct node_gra *net);
gsl_vector *LaplacianSpectrum(struct node_gra *net);
double Synchronizability(struct node_gra *net);

#endif /* !RGRAPH_MISSING_H */
