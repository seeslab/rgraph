/*
  main_tools.c
  $LastChangedDate$
  $Revision$
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "tools.h"

#define EPS 1.e-8

int main()
{
  int result = 0;
  double *vec = allocate_d_vec(4);
  double m, s, mi, ma;

  vec[0] = 1.0;
  vec[1] = 1.0;
  vec[2] = 2.0;
  vec[3] = 4.0;
  
  m = mean(vec, 4);
  s = stddev(vec, 4);
  mi = min(vec, 4);
  ma = max(vec, 4);

  if (fabs(m - 2.0) > EPS)
    result = 1;
  if (fabs(s - sqrt(3.0 / 2.0)) > EPS)
    result = 1;
  if (fabs(mi - 1.0) > EPS)
    result = 1;
  if (fabs(ma - 4.0) > EPS)
    result = 1;

  fprintf(stdout, "mean=%lf stddev=%lf min=%lf max=%lf\n",
	  m, s, mi, ma);

  free_d_vec(vec);

  return result;
}
