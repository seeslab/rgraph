#include <stdlib.h>
#include <math.h>
#include <search.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#include "fillpartitions.h"
#include "io.h"

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
  weightsum *= 2;

  //printf("weightsum: %f\n",weightsum );
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

  // int k;
  // for (i = 0; i < N; i++) {
  //   printf("Node %d (STR:%f ), NEIGH from %d to %d;\n",i,part->nodes[i]->strength,adj->idx[i],adj->idx[i+1]);
  //   for (k=adj->idx[i]; k < adj->idx[i+1]; k++) {
  //     printf("\t %d --> %d (STR:%f) \n",k, adj->neighbors[k], adj->strength[k]);
  //   }
  // }

  free(degree);
  return 0;
}


/**
Compare Edges structures (used with Qsort)

This function compare two edges structures starting by the target node
(if they are equal it compare the source node).

It returns 0 if the edges are identiqual or -1 (resp 1) if the first edge has
smaller (resp. bigger) node ids.
**/
static int
EdgeCompare(const void *p1, const void *p2)
{
        unsigned int y1 = ((Edge*)p1)->node2;
        unsigned int y2 = ((Edge*)p2)->node2;
        if (y1 < y2)
          return -1;
        else if (y1>y2)
          return 1;
        y1 = ((Edge*)p1)->node1;
        y2 = ((Edge*)p2)->node1;
        if (y1 < y2)
          return -1;
        else if (y1>y2)
          return 1;
        return 0;
}



/**
Bipartite projection according to the second column
**/
unsigned int
ProjectBipartEdgeList(unsigned int *nd_in, unsigned int *nd_out, double *weights, int E,
                      Partition **part_p, AdjaArray **adj_p ){
  Partition *part = NULL;
  AdjaArray *adj = NULL;

  // Sort the edges by the nodes to project with.
  int i, count = E;
  unsigned int N = 0;
  Edge *ed = malloc(count*sizeof(Edge));
  for (i = 0; i < count; i++) {
    if (nd_in[i]+1 > N) N = nd_in[i]+1;
    ed[i].node1 = nd_in[i];
    ed[i].node2 = nd_out[i];
    ed[i].strength = weights[i];
  }
  qsort(ed, count, sizeof(Edge), EdgeCompare);

  int Ngroups = N;
  part = CreatePartition(N,Ngroups);

  double fac1 = 0;
  double fac2 = 0;
  double strength;

  // Count the degree of each node and assign a first data structure for the
  // projected network.
  unsigned int degree = 0, Emax = 0;
  for (i = 0; i <= count-1; i++) {
    degree++; // degree of current team
    strength += ed[i].strength; // Strength of current team
    part->nodes[ed[i].node1]->strength += ed[i].strength; // degree of current actor
    if (ed[i].node2 != ed[i+1].node2 || i == count-1) {
      Emax += (degree)*(degree-1) / 2;
      fac2 += strength;
      fac1 += strength * (strength-strength/degree);
      degree = 0;
      strength = 0;
    }
  }

  fac1 = 1. / fac1; // 1/sum_a[m_a(m_a-epsilon)]
  fac2 = 1. / fac2; // 1/(sum_a[m_a])

  Edge *projected = malloc(Emax*sizeof(Edge));
  // Generate all projected edges
  degree = 0;
  unsigned int j = 0, x, y, x0 = 0;
  for (i = 0; i <= count-1; i++) {
    degree++;
    if (ed[i].node2 != ed[i+1].node2 || i == count-1) {
      for (x = x0; x < i+1; x++){
        for (y = x0; y < x; y++) {
          if (ed[x].node1 < ed[y].node1){
            projected[j].node1 = ed[x].node1;
            projected[j].node2 = ed[y].node1;
          }else{
            projected[j].node1 = ed[x].node1;
            projected[j].node2 = ed[y].node1;
          }
          projected[j].strength = ed[x].strength * ed[y].strength;
          j++;
        }
      }
      degree = 0;
      x0 = i+1;
  }
  }
  free(ed);
  // Count the number of unique edges
  qsort(projected, j, sizeof(Edge), EdgeCompare);
  E = j;
  unsigned int *degree_proj = calloc(N,sizeof(unsigned int));
  for (i = 0; i < j; i++) {
    if (EdgeCompare(&projected[i],&projected[i+1]) != 0 || i == j-1) {
      degree_proj[projected[i].node1]++;
      degree_proj[projected[i].node2]++;
    }
    else{
      E --;
  }
  }
  //printf("Create adja array %d nodes / %d edges\n", N,E );

  adj = CreateAdjaArray(N,E);
  unsigned int k = 0;
  unsigned int *idx = malloc(N*sizeof(unsigned int));
  for (i=0; i<N; i++){
    idx[i] = k;
    adj->idx[i] = k;
    k += degree_proj[i];
    //printf("Node %d (degree %d) - idx: %d \n",i,  degree_proj[i], idx[i]);
  }

  //printf("Projected edges:\n");

  for (i = 0; i < j; i++) {
    adj->strength[idx[projected[i].node1]] += projected[i].strength;
    adj->strength[idx[projected[i].node2]] += projected[i].strength;
    if (EdgeCompare(&projected[i],&projected[i+1]) != 0){
      //printf("%d %d %f:\n",projected[i].node1,projected[i].node2,projected[i].strength);
      adj->neighbors[idx[projected[i].node1]] = projected[i].node2;
      adj->neighbors[idx[projected[i].node2]] = projected[i].node1;
      idx[projected[i].node1]++;
      idx[projected[i].node2]++;
    }
  }

  for (i = 0; i < N; i++) {
    //printf("Node %d (STR:%f ), NEIGH from %d to %d;\n",i,part->nodes[i]->strength,adj->idx[i],adj->idx[i+1]);
    part->nodes[i]->strength *= fac2;
    for (k=adj->idx[i]; k < adj->idx[i+1]; k++) {
      adj->strength[k] *= fac1;
      //printf("\t %d --> %d (STR:%f) \n",k, adj->neighbors[k], adj->strength[k]);
    }
  }

  free(degree_proj);
  free(idx);
  free(projected);
  //printf("Setting pointers\n");
  *part_p = part;
  *adj_p = adj;
  return E;
}


/**
Given a partition, dispatch the nodes into modules.
If there are as many modules as nodes, each node goes in one module
exactly. Otherwise, the nodes are dispatched at random between the
modules (without checking if there are empty modules left).
This must be done AFTER initializing the strenght of the nodes,
otherwise, the initial module strength will be incorrect !
**/
void
AssignNodesToModules(Partition *part, gsl_rng *gen){
  unsigned int i,j;
  // if there is as many modules as nodes, assign each node to a module.
  if(part->N == part->M){
	part->nempty = 0;
	for (i=0; i<part->N; i++){
	  part->nodes[i]->module = i;
	  part->modules[i]->size = 1;
	  part->modules[i]->strength = part->nodes[i]->strength;
	  part->modules[i]->first = part->nodes[i];
	  part->modules[i]->last = part->nodes[i];
	}
  }
  // Otherwise dispatch the nodes at random.
  else{
	for (i=0; i<part->N; i++){
	  j = gsl_rng_uniform_int(gen,part->M);
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
	}
  }
}
