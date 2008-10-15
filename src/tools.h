/*
  tools.h
  $LastChangedDate$
  $Revision$
*/

#ifndef RGRAPH_TOOLS_H
#define RGRAPH_TOOLS_H 1

#include "prng.h"

/*
  ---------------------------------------------------------------------
  Vector and matrix memory management
  ---------------------------------------------------------------------
*/
int **allocate_i_mat(int nrows, int ncolumns);
int *allocate_i_vec(int nelem);
double *allocate_d_vec(int nelem);
double **allocate_d_mat(int nrows, int ncolumns);
void free_i_mat(int **data, int nrows);
void free_i_vec(int *data);
void free_d_vec(double *data);
void free_d_mat(double **data, int nrows);

/*
  ---------------------------------------------------------------------
  Random number generation and distributions
  ---------------------------------------------------------------------
*/
int geometric_dist_val(double p, struct prng *gen);

/*
  ---------------------------------------------------------------------
  Statistics
  ---------------------------------------------------------------------
*/
double mean(double *data, int N);
double stddev(double *data, int N);
double max(double *data, int N);
double min(double *data, int N);

/*
  ---------------------------------------------------------------------
  Other
  ---------------------------------------------------------------------
*/
long int fact(long int a);
long double Choose(int a, int b);
long double LogChoose(int a, int b);


#endif /* !RGRAPH_TOOLS_H */
