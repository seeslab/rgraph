/**
@file sannealing.h
@author Guilhem Doulcier
@date 2015
@license GPLv3+
@brief Simulated Annealing (Adjacency Arrays implementation)

This file contains a modularity optimization by simulated annealing
implementation based on adjacency arrays (as defined in
partitions.c). Its main advantage against the old implementation are:

- It is relatively faster (even if it is keeping the same algorithmic
complexity, and thus scaling law). 

- It use a single GeneralSA function to do the four possible SA
  (weighted/unweighted, bipartite or not). The only difference being
  the way the adjacency array is initialized (see fillpartitions.c).
**/
#include <gsl/gsl_rng.h>
#include "partition.h"

unsigned int
GeneralSA(Partition **ppart, AdjaArray *adj,
		  double fac,
		  double Ti, double Tf, double Ts,
		  double cluster_prob,
		  unsigned int nochange_limit,
		  gsl_rng *gen);

unsigned int
SplitModuleSA(unsigned int target, unsigned int empty,
			  double Ti, double Tf, double Ts,
			  unsigned int nochange_limit,
			  Partition *part, AdjaArray *adj,
			  gsl_rng *gen);
