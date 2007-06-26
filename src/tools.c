/*
  tools.c
  $LastChangedDate$
  $Revision$
*/

#include <stdio.h>
#include <stdlib.h>

#include "tools.h"

int **
allocate_i_mat(int nrows, int ncolumns)
{
  int **array;
  int i;

  array = (int**)malloc(nrows * sizeof(int *));
  if(array == NULL){
      fprintf(stderr, "out of memory\n");
      return NULL;
  }

  for(i = 0; i < nrows; i++){
    array[i] = (int*)malloc(ncolumns * sizeof(int));
      if(array[i] == NULL)
	{
	  fprintf(stderr, "out of memory\n");
	  return NULL;
	}
  }

  return array;
}


int *
allocate_i_vec(int nelem)
{
  int *array;

  array = (int*)malloc(nelem * sizeof(int));
  if(array == NULL){
      fprintf(stderr, "out of memory\n");
      return NULL;
  }

  return array;
}


double *
allocate_d_vec(int nelem)
{
  double *array;

  array = (double*)malloc(nelem * sizeof(double));
  if(array == NULL){
      fprintf(stderr, "out of memory\n");
      return NULL;
  }

  return array;
}


double **
allocate_d_mat(int nrows, int ncolumns)
{
  double **array;
  int i;

  array = (double**)malloc(nrows * sizeof(double *));
  if(array == NULL){
      fprintf(stderr, "out of memory\n");
      return NULL;
  }

  for(i = 0; i < nrows; i++){
    array[i] = (double*)malloc(ncolumns * sizeof(double));
      if(array[i] == NULL)
	{
	  fprintf(stderr, "out of memory\n");
	  return NULL;
	}
  }

  return array;
}


void 
free_i_mat(int **data, int nrows)
{
  int i;

  for(i=0; i<nrows; i++){
    free(data[i]);
  }

  free(data);
}

void 
free_i_vec(int *data)
{
  free(data);
}

void
free_d_vec(double *data)
{
  free(data);
}

void
free_d_mat(double **data, int nrows)
{
  int i;

  for(i=0; i<nrows; i++){
    free(data[i]);
  }

  free(data);
}
