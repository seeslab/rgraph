/*
  tools.h
  $LastChangedDate$
  $Revision$
*/

#ifndef RGRAPH_TOOLS_H
#define RGRAPH_TOOLS_H 1

#include <gsl/gsl_rng.h>

/*
  ---------------------------------------------------------------------
  Vector and matrix memory management
  ---------------------------------------------------------------------
*/
int **allocate_i_mat(int nrows, int ncolumns);
void fprintf_i_mat(FILE *outf, int **mat, int nrows, int ncolumns);
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
int geometric_dist_val(double p, gsl_rng *gen);

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
double LogChoose(int a, int b);
double **InitializeFastLogChoose(int LogChooseListSize);
void FreeFastLogChoose(double **LogChooseList,
		       int LogChooseListSize);
double FastLogChoose(int a,
		     int b,
		     double **LogChooseList,
		     int LogChooseListSize);

double *InitializeFastLog(int LogListSize);
void FreeFastLog(double *LogList);
double FastLog(int a, double *LogList, int LogListSize);

double *InitializeFastLogFact(int LogFactListSize);
void FreeFastLogFact(double *LogFactList);
double FastLogFact(int a, double *LogFactList, int LogFactListSize);
double *InitializeFastLogGamma(int LogGammaListSize, double x);
void FreeFastLogGamma(double *LogGammaList);
double FastLogGamma(int r, double *LogGammaList,
		    int LogGammaListSize, double x);

double *InitializeHarmonicList(int HarmonicListSize);
void FreeHarmonicList(double *HarmonicList);

/*
  ---------------------------------------------------------------------
  File operations
  ---------------------------------------------------------------------
*/
int CountLinesInFile(char *inFileName);

#endif /* !RGRAPH_TOOLS_H */
