/**
@file io.h
@author Guilhem Doulcier
@date 2015
@license GPLv3+
@brief Read graph from file and output modules and modularity. 

This file contains all the basic i/o functions used in netcarto. 
**/
#ifndef IO_H__
#define IO_H__
#include <stdlib.h>
#include "partition.h"
#define MAX_LABEL_LENGTH 100

int
EdgeListFileInput(FILE *inFile, int weighted, int bipartite, unsigned int **nodes_in_p,
          unsigned int **nodes_out_p,double **weights_p,
          char ***labels_p, unsigned int *E_p, unsigned int *N_p);
void
TabularOutput(FILE *outf,
       			  char **labels,
       			  Partition *part,
       			  double *connectivity,
       			  double *participation);

void
ClusteringOutput(FILE *outf,
		 Partition *part,
		 char **labels);

int
AssignNodesToModulesFromFile(FILE *inF,
             Partition *part,
             char **labels);

struct node {
   char label[MAX_LABEL_LENGTH];
   unsigned int id;
};
#endif
