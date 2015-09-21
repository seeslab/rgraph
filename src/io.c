#include <stdlib.h>
#include <math.h>
#include <search.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#include "partition.h"
#include "io.h"


/**
  Returns the role number given a (P,z) tuple.

  Cf. Guimera & Amaral, Nature (2005) for the reasoning behind the
  boundaries values.
**/
int GetRole(double P, double z){
  int dest_group;
  if (z < 2.5) {  /* Node is not a hub */
	if (P < 0.050)
	  dest_group = 0;
	else if (P < 0.620)
	  dest_group = 1;
	else if (P < 0.800)
	  dest_group = 2;
	else
	  dest_group = 3;
  }
  else {  /* Node is a hub */
	if (P < 0.300)
	  dest_group = 4;
	else if (P < 0.750)
	  dest_group = 5;
	else
	  dest_group = 6;
  }
  return dest_group;
}


void
TabularOutput(FILE *outf,
			  char **labels,
			  Partition *part,
			  double *connectivity,
			  double *participation){
  unsigned int i=0;
  int rolenb;
  fprintf (outf, "%-30s\tModule\tConnectivity\tParticipation\tRole\n","Label");
  for (i=0;i<part->N;i++){
	rolenb = GetRole(participation[i],connectivity[i])+1;
	fprintf (outf, "%-30s\t%d\t%f\t%f\tR%d\n",
			 labels[i],
			 part->nodes[i]->module,
			 connectivity[i],participation[i],
			 rolenb);
  }
}

void
ClusteringOutput(FILE *outf,
				 Partition *part,
				 char **labels){
  unsigned int mod;
  Node *node;
  for(mod=0;mod<part->M;mod++){
	for(node = part->modules[mod]->first;node!=NULL; node=node->next){
	  fprintf(outf,"%s\t",labels[node->id]);
	}
	fprintf(outf,"\n");
  }
}

int
AssignNodesToModulesFromFile(FILE *inF,
							 Partition *part,
							 char **labels){
  char label[MAX_LABEL_LENGTH];
  char sep[2];
  int j = 0, nfields = 0, nnode = part->N;
  hcreate(part->N);
  int i;
  ENTRY e, *ep;

  for (i=0;i<nnode;i++){
	e.key = labels[i];
	e.data = (void *) i;
	ep = hsearch(e, ENTER);
  }

  while (!feof(inF)){
	nfields=fscanf(inF,"%[^\t\n]%[\t\n]",&label,&sep);
	if (nfields) {
	  e.key = label;
	  ep = hsearch(e, FIND);
	  i = (int) ep->data;
	  nnode--;

	  if (!part->modules[j]->size){
		part->nempty --;
		part->nodes[i]->module = j;
		part->modules[j]->size = 1;
		part->modules[j]->strength = part->nodes[i]->strength;
		part->modules[j]->first = part->nodes[i];
		part->modules[j]->last = part->nodes[i];
	  }
	  else{
		part->nodes[i]->module = j;
		part->modules[j]->size++;
		part->modules[j]->strength += part->nodes[i]->strength;
		part->modules[j]->last->next = part->nodes[i];
		part->nodes[i]->prev = part->modules[j]->last;
		part->modules[j]->last = part->nodes[i];
	  }
	  if(sep[0]=='\n')
		j++;
	}
  }

  CompressPartition(part);
  hdestroy();
  return nnode;
}

/**
Normalize edges weight and node strength and store them in the
Partition and AdjaArray structures.

This function assume that the edge list is UNDIRECTED and WITHOUT
duplicates.

If W_ij is the adjacency matrix and W=\sum_i\sum_j W_ij the sum of its
elements. The normalized strength is:
- For an edge:  A_ij = W_ij / W
- For a node: k_i = sum_i (W_ij) / W

@param nd_in,nd_out,weight Edge list, nodes id and weight of each edge.
@param adj The adjacency array to fill.
@param part The partition to fill.
@param normalize If set to 0, the strength and weight are not normalized by W.
**/
int
EdgeListToAdjaArray(int *nd_in, int *nd_out, double *weight,
					AdjaArray *adj, Partition *part, int normalize){

  unsigned int N = adj->N, E = adj->E;
  double weightsum = 0;
  unsigned int i;

  int *degree = (int*) calloc(N,sizeof(int));
  if (degree == NULL){
	   return 1;
  }

  // Compute degrees and the sum of edges weight.
  for (i=0;i<E;i++){
  	weightsum += weight[i];

  	part->nodes[nd_in[i]]->strength +=  weight[i];
  	part->nodes[nd_out[i]]->strength +=  weight[i];

  	degree[nd_in[i]] += 1;
  	degree[nd_out[i]] += 1;
  }

  // If the weights and strength are not to be normalized, set the
  // normalization constant to 1 (neutral element for division).
  if (!normalize)
	 weightsum = 1;

  // Set the start of the index of the neighbors and store the
  // normalized strength.
  int pos = 0;
  for (i=0;i<N;i++){
  	adj->idx[i] = pos;
  	pos += degree[i];
  	part->nodes[i]->strength /=  weightsum;
  }

  // Fill the edges properties (target and normalized weights).  We
  // re-use the degree array to keep track of the position to fill for
  // each node. The first edge encountered for a node will be stored
  // in the first position for this node:
  // (idx[nodeid+1]-degree[node]=idx[nodeid]). Successives edges will
  // be stored in successive positions (because degree[node] is
  // decremented each time).
  int idx_in, idx_out;
  for (i=0;i<E;i++){
  	// Get the positions to fill.
  	idx_in  = adj->idx[nd_in[i]+1]-degree[nd_in[i]];
  	idx_out = adj->idx[nd_out[i]+1]-degree[nd_out[i]];

  	// Fill them with target and strength.
  	adj->neighbors[idx_in] = nd_out[i];
    adj->strength[idx_in] = weight[i]/weightsum;
  	adj->neighbors[idx_out] = nd_in[i];
  	adj->strength[idx_out] = weight[i]/weightsum;

  	// Update the counter to fill the next position correctly.
  	degree[nd_in[i]]--;
  	degree[nd_out[i]]--;

  }
  return 0;
}
