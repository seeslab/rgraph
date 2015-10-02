#include "partition.h"
#include "math.h"
#include <stdio.h>

/**
Allocate the memory needed for an adjacency array
@param N number of nodes
@param E number of edges
 **/
AdjaArray *
CreateAdjaArray(unsigned int N, unsigned int E){
  AdjaArray *adj;
  adj = (AdjaArray *) malloc(sizeof(AdjaArray));
  if (adj==NULL)
	perror("Error while allocating adjacency array");
  adj->N = N;
  adj->E = E;
  adj->idx = (unsigned int *) calloc(N+1,sizeof(unsigned int ));
  adj->neighbors = (unsigned int *) calloc(2*E+1,sizeof(unsigned int ));
  adj->idx[N] = 2*E;
  adj->strength = (double *) calloc(2*E+1,sizeof(double));
  if (adj->idx==NULL || adj->neighbors == NULL || adj->strength == NULL)
	perror("Error while allocating adjacency array elements");
  return adj;
}

/**
Free the memory used by an adjacency array.
**/
void
FreeAdjaArray(AdjaArray *adj){
  free(adj->idx);
  free(adj->neighbors);
  free(adj->strength);
  free(adj);
  adj = NULL;
}

/**
Allocate the memory for a partition (nodes included).
@param N number of nodes
@param M number of modules
 **/
Partition *
CreatePartition(unsigned int N, unsigned int M){
  Partition *part = NULL;
  int i;
  part = malloc(sizeof(Partition));
  if (part==NULL)
	perror("Error while allocating partition");

  part->N = N;
  part->M = M;
  part->nempty = M;
  part->nodes = (Node **) malloc(N*sizeof(Node*));
  part->modules = (Module **) malloc(N*sizeof(Module*));
  if (part->nodes==NULL || part->modules == NULL)
	perror("Error while allocating partition component");


  for(i=0;i<N;i++){
	part->nodes[i] = malloc(sizeof(Node));
	if (part->nodes[i]==NULL)
	  perror("Error while allocating node");

	part->nodes[i]->id = i;
	part->nodes[i]->module = 0;
	part->nodes[i]->strength = 0;
	part->nodes[i]->next = NULL;
	part->nodes[i]->prev = NULL;
  }
  for(i=0;i<M;i++){
	part->modules[i] = malloc(sizeof(Module));
	if (part->modules[i]==NULL)
	  perror("Error while allocating module");

	part->modules[i]->id=i;
	part->modules[i]->strength = 0;
	part->modules[i]->size=0;
	part->modules[i]->first = NULL;
	part->modules[i]->last = NULL;
  }
  return(part);
}

/**
 Deep copy of a partiton
**/
Partition *
CopyPartitionStruct(Partition *part){
  Partition *copy = NULL;
  int i;
  copy = CreatePartition(part->N,part->M);
  copy->nempty = part->nempty;

  // Copy the nodes...
  for(i=0;i<part->N;i++){
	copy->nodes[i]->id = part->nodes[i]->id;
	copy->nodes[i]->module = part->nodes[i]->module;
	copy->nodes[i]->strength = part->nodes[i]->strength;
  }

  // Copy the modules...
  for(i=0;i<part->M;i++){
	copy->modules[i]->id = part->modules[i]->id;
	copy->modules[i]->strength = part->modules[i]->strength;
	copy->modules[i]->size = part->modules[i]->size;
  }

  // Copy the links in the doubly linked list of nodes (=modules).
  for(i=0;i<part->N;i++){
	if(part->nodes[i]->next != NULL)
	  copy->nodes[i]->next = copy->nodes[part->nodes[i]->next->id];
	else
	  copy->nodes[i]->next = NULL;
	if(part->nodes[i]->prev != NULL)
	  copy->nodes[i]->prev = copy->nodes[part->nodes[i]->prev->id];
	else
	  copy->nodes[i]->prev = NULL;
  }
  for(i=0;i<part->M;i++){
	if (part->modules[i]->first != NULL)
	  copy->modules[i]->first = copy->nodes[part->modules[i]->first->id];
	else
	  copy->modules[i]->first = NULL;
	if (part->modules[i]->last != NULL)
	  copy->modules[i]->last = copy->nodes[part->modules[i]->last->id];
	else
	  copy->modules[i]->last = NULL;
  }
  return(copy);
}

/**
Free the memory used by a partition
**/
void
FreePartition(Partition *part){
  int i;
  for (i=0;i<part->N;i++)
	free(part->nodes[i]);
  for (i=0;i<part->M;i++)
	free(part->modules[i]);
  free(part->nodes);
  free(part->modules);
  free(part);
  part = NULL;
}


/**
Compute the modularity roles metrics of all nodes.

@param part the Partition
@param adj the adjacency array
@param connectivity (must point to an array of size N). The within module degree z-score.
@param participation (must point to an array of size N). The evenness of module degrees.
 **/
void
PartitionRolesMetrics(Partition *part, AdjaArray *adj, double *connectivity, double *participation){
  unsigned int N,M,i,j, mod;
  double *strengthToModule,*mean,*std;

  // Number of nodes and modules.
  N = adj->N;
  M = part->M;

  ////Compute all strengths from a node to a module.
  strengthToModule = (double*) calloc(N*M,sizeof(double));
  if (strengthToModule==NULL)
	  perror("Error while computing roles metrics");

  for (i=0; i<N; i++){
	for (j=adj->idx[i]; j<adj->idx[i+1]; j++){
	  mod = part->nodes[adj->neighbors[j]]->module;
	  strengthToModule[i+N*mod] += adj->strength[j];
	}
  }

  //// Connectivity is the within module z-score
  //// C_i = X-E(X)/STD(X) with X = K_iG if i is in module G.
  mean = (double*) calloc(M,sizeof(double));
  std = (double*) calloc(M,sizeof(double));
  if (mean==NULL || std ==NULL)
	  perror("Error while computing roles metrics");

  // We compute the mean µ and then the std(X-µ)=std(X) to acheive
  // numerical stability.
  for (i=0; i<N; i++){
	mod = part->nodes[i]->module;
	mean[mod] += strengthToModule[i+N*mod];
  }
  for (j=0; j<M; j++)
	mean[j] /= part->modules[j]->size;

  for (i=0; i<N; i++){
	mod = part->nodes[i]->module;
	std[mod] += (strengthToModule[i+N*mod]-mean[mod]) * (strengthToModule[i+N*mod]-mean[mod]);
  }
  for (j=0; j<M; j++)
	std[j] = sqrt(std[j]/part->modules[j]->size);

  // Compute the z-score.
  for (i=0; i<N; i++){
	mod = part->nodes[i]->module;
	if (std[mod])
	  connectivity[i] = (strengthToModule[i+N*mod] - mean[mod]) / std[mod];
	else
	  connectivity[i] = 0;
  }


  // Participation is the 1-simpson index of node to modules strength.
  // P_i = 1 - SUM(K_iG^2/k_i^2).
  for (i=0; i<N; i++){
	// Participation coefficient of unconnected nodes is null.
	if(adj->idx[i]==adj->idx[i+1]) //Node i is of null degree.
	  participation[i] = 0;
	else{
	  //we cannot use part->nodes[i]->strength for denominator because
	  //it is not correct for bipartite network thanks to the
	  //-epsilon.
	  double strength = 0;
	  for (j=0;j<M;j++){
		strength += strengthToModule[i+N*j];
		participation[i] += strengthToModule[i+N*j]*strengthToModule[i+N*j];
	  }
	  participation[i] = 1.0 - participation[i]/(strength*strength);
	}
  }

  //free memory.
  free(mean);
  free(std);
  free(strengthToModule);
}

/**
Compute the modularity of a partition.

The modularity is given by:
M(P) = \sum[i!=j,P(i)==P(j)] A_ij - k_ik_j

If diagonal term evaluate to true, the formula is:
M(P) = \sum[(i,j),P(i)==P(j)] A_ij - k_ik_j

**/
double
PartitionModularity(Partition *part, AdjaArray *adj, int diagonal_term){
  int mod;
  Node *node1, *node2;
  double modularity = 0.0;
  double aij;
  int i;

  // For all modules...
  for(mod=0;mod<part->M;mod++){

	// For all pair of node i!=j in a module...
	for(node1 = part->modules[mod]->first;node1!=NULL; node1=node1->next){
	  for(node2 = node1->next ; node2!=NULL; node2=node2->next){

		// Get the value of the edge (if they are connected)
		aij = 0;
		for (i=adj->idx[node1->id]; i<adj->idx[node1->id+1]; i++){
		  if (adj->neighbors[i] == node2->id){
			aij = adj->strength[i];
			break; // found the edge, we can stop.
		  }
		}

		// Update modularity.
		modularity += 2 * (aij - node1->strength*node2->strength);
	  }
	}
  } // End looping through modules.

  // If necessary remove the diagonal term (only switch the modularity
  // by a constant as a node is always in the same module as itself).
  if (diagonal_term){
	for(i=0 ; i<part->N; i++)
	  modularity -= part->nodes[i]->strength*part->nodes[i]->strength;
  }

  return(modularity);
}

/**
Remove empty modules from a partition. Consolidating their indices so
that all modules have consecutive indices in 1...M.
**/
void
CompressPartition(Partition *part){
  Module **newmodules;
  unsigned int M, i,j=0, *empty_id;
  Node *node;

  // If there is no empty modules, do nothing.
  if (!part->nempty) return;

  // The new number of modules is M.
  M = part->M - part->nempty;
  newmodules =  (Module **) malloc(M*sizeof(Module*));
  if (newmodules==NULL)
	  perror("Error while compressing partition");


  // Free the empty modules and store their ids.
  empty_id = (unsigned int *) calloc(part->nempty,sizeof(unsigned int));
  if (empty_id==NULL)
	perror("Error while compressing partition");

  for (i=0;i<part->M;i++){
	if (!part->modules[i]->size){
	  empty_id[j] = i;
	  j++;
	  free(part->modules[i]);
	  part->modules[i]=NULL;
	}
  }

  // Starting from the end, loop through the modules and move the ones
  // that were not freed (!=NULL) to one of the empty ids in the
  // begining of the list.
  j = 0;
  for (i=(part->M)-1;i>=M;i--){
	if (part->modules[i]!=NULL) {
	  for(node=part->modules[i]->first; node!=NULL; node = node->next){
		node->module = empty_id[j];
	  }
	  part->modules[empty_id[j]] = part->modules[i];
	  part->modules[empty_id[j]]->id = empty_id[j];
	  j++;
	}
  }
  free(empty_id);

  // Shorten the module index stored by the partition struct.
  for (i=0;i<M;i++)
	newmodules[i] = part->modules[i];
  free(part->modules);
  part->modules = newmodules;

  // Settting the Partiton properties.
  part->nempty = 0;
  part->M = M;
}

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
