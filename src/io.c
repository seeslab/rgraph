#include <stdlib.h>
#include <math.h>
#include <search.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#include "partition.h"

typedef struct EdgeList {
  unsigned int N; //! Number of nodes.
  unsigned int E; //! Number of edges.
  unsigned int *node_in; //!< 
  unsigned int *node_out; //!<
  char **labels; //!< node labels
  double *weight; //!< weight of the corresponding edge
} EdgeList;

#define MAX_LABEL_LENGTH 100
EdgeList* CreateEdgeList(int N, int E, int labels); // Allocate the memeory for an edge list.
void FreeEdgeList(EdgeList *elist); // Free an edge list.
EdgeList* CopyEdgeList(EdgeList *elist); // Deep copy an edge list
EdgeList* EdgeListFromFile(FILE *inF); // Read a file and build the edge list
int MCMCShuffleEdgeList(EdgeList elist,int bipartite); // Shuffle edges using a MCMC heuristic.  

// Project a bipartite network.
int ProjectBipartite(AdjaArray* adj, Partition **projected_part,
					 AdjaArray **projected_adj);


// Convert an EdgeList into an AdjaArray and Partition
int EdgeListToAdjaArray(int *nd_in, int *nd_out, double *weight,
						AdjaArray *adj, Partition *part, int normalize);
// Read a partition from a file. 
int AssignNodesToModulesFromFile(FILE *inF, Partition *part, char **labels);

// Output functions
void TabularOutput(FILE *outf, char **labels, Partition *part,
				   double *connectivity, double *participation);
void ClusteringOutput(FILE *outf,  Partition *part, char **labels);



EdgeList*
CreateEdgeList(int N, int E, int labels){
  EdgeList* elist=NULL;
  if (elist==NULL){
	perror("Error while allocating edgelist");
	return NULL;
  }
  elist->N = N;
  elist->E = E;
  elist->node_in  = (unsigned int *) calloc(E,sizeof(unsigned int));
  elist->node_out = (unsigned int *) calloc(E,sizeof(unsigned int));
  elist->weight = (double *) calloc(E,sizeof(double));

  if (elist->node_in==NULL || elist->neighbors == NULL || elist->weight == NULL){
	perror("Error while allocating edge list elements");
	return NULL;
  }

  if (!labels)
	elist->labels = NULL;
  else{
	elist->labels = (char**) malloc(N*sizeof(char*));
	
	if (elist->labels==NULL){
		perror("Error while allocating EdgeList label list");
		return NULL;
	} 

	for (i=0;i<N;i++){
	  elist->labels[i] = (char *) malloc(MAX_LABEL_LENGTH * sizeof(char));
	  if (elist->labels[i]==NULL){
		perror("Error while allocating EdgeList label list elements");
		return NULL;
	  }
	}
  }
  return elist;
}

void
FreeEdgeList(EdgeList *elist){
  free(elist->node_in);
  free(elist->node_out);
  free(elist->weight);
  if (elist->lablels != NULL)
	free(elist->labels);
  free(elist);
  elist = NULL;
}

EdgeList*
CopyEdgeList(EdgeList elist){
  int labels = 0;
  int i;
  if (elist->labels != NULL)
	labels = 1;
  
  EdgeList *copy = CreateEdgeList(elist->N,elist->E,labels);
  for (i=0;i<elist->E;i++){
	copy->node_in[i] = elist->node_in[i];
	copy->node_out[i] = elist->node_out[i];
	copy->weight[i] = elist->weight[i];
  }

  if (elist->labels != NULL){
	for (i=0;i<elist->N;i++)
	  strcpy(copy->labels[i],elist->labels[i]);			 
  }
}

EdgeList*
EdgeListFromFile(FILE *inF, int weighted)
{
  char *line = NULL;
  size_t bufsiz = 0;
  ssize_t nbytes;
  char label1[MAX_LABEL_LENGTH], label2[MAX_LABEL_LENGTH];
  double weight;
  
  // Go through the input file
  while ((nbytes = getline(&line, &bufsiz, inF)) != -1){
    /* Read the labels (and weight, if necessary) */
    if (weighted == 0) {
      noReadItems = sscanf(line, "%s %s", label1, label2);
	  if(noReadItems != 2){
		printf ("Failed to read input: not enough fields in line %s (%d!=2). \n",line, noReadItems);
		return NULL;
	  }
	  weight = 1;
    }
    else {
      noReadItems = sscanf(line,"%s %s %lf", label1, label2, &weight);
	  if(noReadItems != 3){
		printf ("Failed to read input: not enough fields in line %s (%d!=3). \n",line, noReadItems);
		return NULL;
	  }
    }

	// TODO ! 
  } // Done Reading
  free(line);
  return elist;
}


/**
Shuffle 


 **/
int
MCMCShuffleAdjaArray(EdgeList elist,int bipartite, int weighted,int success); 

/**
@param adj
@param projected_part
@param projected_adj
@parm first_team_id 
 **/
int
ProjectBipartite(AdjaArray* adj, Partition **projected_part,
				 AdjaArray **projected_adj,unsigned int first_team_id,
				 unsigned int side, unsigned int M)
{
  unsigned int N, E,start,stop;
  if (side){
	N = first_team_id;
	start = first_team_id;
	stop = adj->N;
  }else{
	N = adj->N - first_team_id;
	start = 0;
	stop = first_team_id;
  }
  
  // Loop through pairs of teams neighbors.
  unsigned int i,j,k,node1,node2;
  for (i = start; i++; i<stop){
	for (j=adj->idx[i]; j<adj->idx[i+1]; j++){
	  node1 = adj->neighbors[j];
	  for (k=j; k<adj->idx[i+1]; k++){
		node2 = adj->neighbors[k];
		// add node1,node2,j
	  }
	}
  }
	
}




void
TabularOutput(FILE *outf,
			  char **labels,
			  Partition *part,
			  double *connectivity,
			  double *participation){
  int i=0;
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
  int mod;
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
  int mod;
  Node *node;
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
	  i = ep->data;
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

/*
, int bipartite
  @param bipartite If true, the nodes in nd_out are considered different
than nodes in nd_in even if they have the same id. 

  if (bipartite){
	unsigned int id_first_team = nd_in[0];
 	for (i=1;i<E;i++){
	  if (nd_in[i]>id_first_team)
		id_first_team = nd_in[i];
	}
	for (i=1;i<E;i++){
	  nd_out[i] += id_first_team; 
  }
*/


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

  int N = adj->N, E = adj->E;
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
	adj->neighbors[idx_out] = nd_in[i];
	adj->strength[idx_out] = weight[i]/weightsum;
	adj->strength[idx_in] = weight[i]/weightsum;

	// Update the counter to fill the next position correctly.
	degree[nd_in[i]]--;
	degree[nd_out[i]]--;
  }
  return 0;
}
