/* $Revision$ */

#ifndef RGRAPH_TOOLS_H
#define RGRAPH_TOOLS_H 1

int **allocate_i_mat(int nrows, int ncolumns);
int *allocate_i_vec(int nelem);
double *allocate_d_vec(int nelem);
double **allocate_d_mat(int nrows, int ncolumns);
void free_i_mat(int **data, int nrows);
void free_i_vec(int *data);
void free_d_vec(double *data);
void free_d_mat(double **data, int nrows);

#endif /* !RGRAPH_TOOLS_H */
