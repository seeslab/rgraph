/*
  matrix.c
  $LastChangedDate: 2008-10-13 19:13:23 -0500 (Mon, 13 Oct 2008) $
  $Revision: 130 $
*/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

#include "tools.h"
#include "graph.h"
#include "matrix.h"

/*
  Build the adjacency matrix of a network
*/
int **
AdjacencyMat(struct node_gra *net)
{
  int nnod=CountNodes(net);
  struct node_gra *p=net;
  struct node_lis *n;
  int **A, i, j;

  /* Create and initialize the matrix */
  A = allocate_i_mat(nnod, nnod);
  for (i=0; i<nnod; i++)
    for (j=0; j<nnod; j++)
      A[i][j] = 0;

  /* Add the ones */
  while ((p = p->next) != NULL) {
    n = p->neig;
    while ((n = n->next) != NULL) {
      A[p->num][n->ref->num] = 1;
    }
  }

  /* Done */
  return A;
}

/*
  Build the Laplacian matrix of a network
*/
int **
LaplacianMat(struct node_gra *net)
{
  struct node_gra *p=net;
  int nnod=CountNodes(net);
  int **L, i, j;;

  /* Get the adjacency matrix */
  L = AdjacencyMat(net);

  /* Make it the Laplacian */
  for (i=0; i<nnod; i++) {
    p = p->next;
    if (L[i][i] > 0) {
      fprintf(stderr, "ERROR: Graph contains self-loops (this will crash)");
      return NULL;
    }
    else {
      L[i][i] = NodeDegree(p);
    }
    for (j=i+1; j<nnod; j++) {
      L[i][j] = -L[i][j];
      L[j][i] = -L[j][i];
    }
  }
  
  /* Done */
  return L;
}

/*
  Get the spectrum of the Laplacian
*/
gsl_vector *
LaplacianSpectrum(struct node_gra *net)
{
  int nnod=CountNodes(net);
  int **L=LaplacianMat(net);
  int i, j;
  gsl_matrix *gslL;
  gsl_eigen_symm_workspace *workspace;
  gsl_vector *eval;

  /* Create a Laplacian ready for GSL use */
  gslL = gsl_matrix_calloc(nnod, nnod);
  for (i=0; i<nnod; i++)
    for (j=0; j<nnod; j++)
      gsl_matrix_set(gslL, i, j, L[i][j]);
  
  /* Get the eigenvalues */
  workspace = gsl_eigen_symm_alloc(nnod);
  eval = gsl_vector_alloc(nnod);
  gsl_eigen_symm(gslL, eval, workspace);
  gsl_sort_vector(eval);

  /* Free memory and done */
  free_i_mat(L, nnod);
  gsl_matrix_free(gslL);
  gsl_eigen_symm_free(workspace);

  return eval;
}

/*
  Get the "synchronizability", that is the ratio between the largest
  and the smallest (non-zero) eigenvalues of the Laplacian.
*/
double
Synchronizability(struct node_gra *net)
{
  gsl_vector *spec=LaplacianSpectrum(net);
  int count=0, nnod=CountNodes(net);
  double l2=-1.0, lN=-1.0;
  
  lN = gsl_vector_get(spec, nnod - 1);
  while (l2 < 1e-10)
    l2 = gsl_vector_get(spec, count++);
  
  gsl_vector_free(spec);
  return lN / l2;
}
