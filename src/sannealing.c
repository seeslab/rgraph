#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "bipartite.h"
#include "sannealing.h"

#define EPSILON_MOD 1.e-6

/**
Compute a the number of iterations for the SA.

@param nnod Number of nodes
@param fac Iteration factor
@param individual_movements A pointer to where to write the number of local movements.
@param collective_movements A pointer to where to write the number of collectives movements.
 **/
void iterationNumber(int nnod,
					 double fac,
					 int *individual_movements,
					 int *collective_movements){
  if (fac * (double)(nnod * nnod) < 10)
	*individual_movements = 10;
  else
	*individual_movements = floor(fac * (double)(nnod * nnod));

  if (fac * (double)nnod < 2)
	*collective_movements = 2;
  else
	*collective_movements = floor(fac * (double)nnod);
}

struct group *
SACommunityIdent(struct node_gra *net,
				 double Ti, double Tf, double Ts,
				 double fac,
				 int ngroup,
				 char initial_sw,
				 int collective_sw,
				 char output_sw,
				 gsl_rng *gen){
  struct group *ret;
  double **cmat;
  unsigned int N;
  N = CountNodes(net);
  cmat = CostMatrix(net);
  /* ret = GeneralSA(cmat, N, net, */
  /* 				  fac, Ti, Tf,  Ts, */
  /* 				  gen); */
  free(cmat);
  return(ret);
}

struct group *
SACommunityIdentWeight(struct node_gra *net,
						 double Ti, double Tf, double Ts,
						 double fac,
						 int merge,
						 gsl_rng *gen)
{
  struct group *ret;
  double **cmat;
  unsigned int N;
  N = CountNodes(net);
  cmat = CostMatrixWeighted(net);
  /* ret = GeneralSA(cmat, N, net, */
  /* 				  fac, Ti, Tf,  Ts, */
  /* 				  gen); */
  free(cmat);
  return(ret);
}

struct group *
SACommunityIdentBipart(struct binet *binet,
						 double Ti, double Tf, double Ts,
						 double fac,
						 int ngroup,
						 char initial_sw,
						 int collective_sw,
						 char output_sw,
						 gsl_rng *gen)
{
  struct group *ret;
  double **cmat;
  struct node_gra *projected ;
  unsigned int N;
  N = CountNodes(binet->net1);
  cmat = CostMatrixBipart(binet);
  projected = ProjectBipart(binet);
  /* ret = GeneralSA(cmat, N, projected, */
  /* 				  fac, Ti, Tf,  Ts, */
  /* 				  gen); */
  free(cmat);
  return(ret);
}

struct group *
SACommunityIdentBipartWeighted(struct binet *binet,
								 double Ti, double Tf, double Ts,
								 double fac,
								 int ngroup,
								 char initial_sw,
								 int collective_sw,
								 char output_sw,
								 gsl_rng *gen)

{
  struct group *ret;
  double **cmat;
  struct node_gra *projected ;
  unsigned int N;
  N = CountNodes(binet->net1);
  cmat = CostMatrixBipartWeighted(binet);
  projected = ProjectBipart(binet);
  /* ret = GeneralSA(cmat, N, projected, */
  /* 				  fac, Ti, Tf,  Ts, */
  /* 				  gen); */
  free(cmat);
  return(ret);
}


struct group*
convertPartitionToGroup(struct partition *part, struct node_gra *net)
{
  struct group *group = NULL, *lgroup = NULL;
  struct node_gra *node;
  struct group **glist;
  unsigned int i;

  // Create the header
  group = CreateHeaderGroup();
  lgroup = group;
  
  // Create a grouplist for faster access.
  glist = (struct group **) calloc(part->nmod, sizeof(struct group *));

  // Create the groups
  for (i=0;i<part->nmod;i++){
	if (part->group[i] != NULL){
	  glist[i] = CreateGroup(lgroup, i);
	  lgroup = lgroup->next;	  
	}
  }

  node = net;
  // Add the nodes to their groups.
  while ((node = node->next) != NULL) {
	i = part->module[node->num];
	AddNodeToGroup(glist[i], node);
  }

  // Free memory
  free(glist);
	
  return(CompressPart(group));
}

double **CostMatrix(struct node_gra *net){
  struct node_gra *p, *q;
  struct node_lis *neigh;
  int nnod = CountNodes(net);
  double **cmat = allocate_d_mat(nnod, nnod);
  double links = 0, fac1=0, fac2=0, s1,s2;
  
  p = net;
  while ((p = p->next) != NULL) {
	links += (double)NodeDegree(p);
  }
  fac1 = 2/links;
  fac2 = 1/(2*links*links);

  p = net;
  while ((p = p->next) != NULL) {
	s1 = (double)NodeDegree(p);
	q = net;
	while ((q = q->next) != NULL) {
	  s2 = (double)NodeDegree(q);
	  if(cmat[q->num][p->num] != 0)
		cmat[p->num][q->num] = cmat[q->num][p->num];
	  else{
		cmat[p->num][q->num] = - fac2 * s1 * s2;
		neigh = p->neig;
		while ((neigh = neigh->next) != NULL){
		  if (neigh->node == q->num){
			cmat[p->num][q->num] += fac1;
			break;
		  }
		}
	  }
	}
  }
  return cmat;
}

double **CostMatrixWeighted(struct node_gra *net){
  struct node_gra *p, *q;
  struct node_lis *neigh;
  int nnod = CountNodes(net);
  double **cmat = allocate_d_mat(nnod, nnod);
  double links = 0, fac1=0, fac2=0, s1,s2;
  
  p = net;
  while ((p = p->next) != NULL) {
	links += (double)NodeStrength(p);
  }
  fac1 = 2/links;
  fac2 = 1/(2*links*links);

  p = net;
  while ((p = p->next) != NULL) {
	s1 = (double)NodeStrength(p);
	q = net;
	while ((q = q->next) != NULL) {
	  s2 = (double)NodeStrength(q);
	  if(cmat[q->num][p->num] != 0)
		cmat[p->num][q->num] = cmat[q->num][p->num];
	  else{
		cmat[p->num][q->num] = - fac2 * s1 * s2;
		neigh = p->neig;
		while ((neigh = neigh->next) != NULL){
		  if (neigh->node == q->num){
			cmat[p->num][q->num] += fac1*neigh->weight;
			break;
		  }
		}
	  }
	}
  }
  return cmat;
}

	
double **CostMatrixBipart(struct binet *binet){
  struct node_gra *p, *q;
  int nnod = CountNodes(binet->net1);
  double **cmat = allocate_d_mat(nnod, nnod);
  double links=0, fac1=0, fac2=0,s1=0,s2=0;
  p = binet->net2;
  
  while ((p = p->next) != NULL) {
	links = (double)NodeDegree(p);
	fac1 += links * (links - 1);
	fac2 += links;
  }
  fac1 = 1. / fac1; // 1/sum_a[m_a(m_a-1)]
  fac2 = 1. / (fac2*fac2); // 1/(sum_a[m_a])**2

  p = binet->net1;
  while ((p = p->next) != NULL) {
	s1 = (double)NodeDegree(p);
	q = binet->net1;
	while ((q = q->next) != NULL) {
	  s2 = (double)NodeDegree(q);
	  if(cmat[q->num][p->num] != 0)
		cmat[p->num][q->num] = cmat[q->num][p->num];
	  else
		cmat[p->num][q->num] = fac1 * (double)NCommonLinksBipart(p, q) - fac2 * s1 * s2;
	}
  }
  return cmat;
}


double **CostMatrixBipartWeighted(struct binet *binet){
  struct node_gra *p, *q;
  int nnod = CountNodes(binet->net1);
  struct node_lis *neigh;
  double **cmat = allocate_d_mat(nnod, nnod);
  double links, fac1, fac2, s1, s2, min;
  p = binet->net2;
  
  while ((p = p->next) != NULL) {
	links = NodeStrength(p);
	min = links;
	neigh = p->neig;
	while ((neigh = neigh->next) != NULL){
	  if (neigh->weight < min){
		min = neigh->weight;
	  }
	}
	fac1 += links * (links - min);
	fac2 += links;
  }
  fac1 = 1. / fac1; // 1/sum_a[m_a(m_a-1)]
  fac2 = 1. / (fac2*fac2); // 1/(sum_a[m_a])**2

  p = binet->net1;
  while ((p = p->next) != NULL) {
	s1 = (double)NodeStrength(p);
	q = binet->net1;
	while ((q = q->next) != NULL) {
	  s2 = (double)NodeStrength(q);
	  if(cmat[q->num][p->num] != 0)
		cmat[p->num][q->num] = cmat[q->num][p->num];
	  else
		cmat[p->num][q->num] = fac1 * (double)SumProductsOfCommonWeightsBipart(p, q) - fac2 * s1 * s2;
	}
  }
  return cmat;
}


struct partition * InitializePartition(unsigned int N){
  struct partition *part = NULL;
  struct llist *element;
  unsigned int i;

  part = (struct partition *) malloc(sizeof(struct partition));
  part->nnod = N;
  part->nmod = N;

  part->module = (unsigned int *) malloc(N*sizeof(unsigned int ));
  part->size = (unsigned int *) malloc(N*sizeof(unsigned int ));
  part->group = (struct llist **) calloc(N,sizeof(struct llist *));
  
  // Assign each node to a module.
  for(i=0;i<N;i++){
	element = (struct llist *) malloc(sizeof(struct llist));
	element->value = i;
	element->next = NULL;
	part->group[i] = element;
	part->module[i] = i;
	part->size[i] = 1;
  }
  return(part);
}

struct partition *CopySimplePartition(struct partition *part){
  struct partition *copy = NULL;
  struct llist *element, *previous;
  unsigned int i, j;
  copy = (struct partition *) malloc(sizeof(struct partition));
  copy->nnod = part->nnod;
  copy->nmod = part->nmod;
  
  copy->module = (unsigned int *) malloc(part->nnod*sizeof(unsigned int ));
  copy->size = (unsigned int *) malloc(part->nnod*sizeof(unsigned int ));
  copy->group = (struct llist **) calloc(part->nmod,sizeof(struct llist *));
  
  for(i=0;i<part->nnod;i++){
	j = part->module[i];
	previous = part->group[j];
	
	// Add the element at the begining of the list.
	element = (struct llist *) malloc(sizeof(struct llist));
	element->value = i;
	element->next = previous;
	copy->group[j] = element;

	// Copy the module number and update the size. 
	copy->module[i] = j;
	copy->size[j] += 1;
  }
  return(copy);
}

void FreeSimplePartition(struct partition *part){
  unsigned int i;
  free(part->module);
  free(part->size);
  for(i=0;i<part->nmod;i++){
	free(part->group[i]);
  }
  free(part);
}

struct group *
GeneralSA(double **cmat, unsigned int N,
		  struct node_gra *net,
		  double fac,
		  double Ti, double Tf, double Ts,
		  gsl_rng *gen)
{
  struct partition *part = NULL, *best_part = NULL;
  struct group *grp;
  unsigned int individual_movements, collective_movements;
  unsigned int target, oldg, newg, i, g1, g2;
  int empty=0;
  double T = Ti;
  unsigned int ngroup = N, nochange_count=0, nochange_limit = 25;
  double dE=0.0, E=0.0, previousE=0.0, best_E=-1.0/0.0; //initial best is -infinity.
  struct llist *element, *element2;
  part = InitializePartition(N);
  
  // Get the number of individual and collective movements. 
  iterationNumber(N, fac, &individual_movements, &collective_movements);

  /// SIMULATED ANNEALING ///
  printf ("#T\tE\n");
  for (T=Ti; T > Tf; T = T*Ts) {
	printf ("%1.5f\t%1.5f\n",T,E);
	//// INDIVIDUAL MOVEMENTS. ////
   	for (i=individual_movements; i; i--) {
	  //// Select a node and a target group.  
	  target = floor(gsl_rng_uniform(gen) * (double)N);
	  oldg = part->module[target];
	  do {
		newg = floor(gsl_rng_uniform(gen) * ngroup);
	  } while (newg == oldg);

	  //// Computing the difference in energy.
	  // Removing the node from its group.
	  dE = 0.0;
	  element = part->group[oldg];
	  while ((element = element->next) != NULL){
		dE -= cmat[target][element->value];
	  }
	  // Adding it to its new group.
	  element = part->group[newg];
	  while ((element = element->next) != NULL) 
		dE += cmat[target][element->value];

	  //// Accept or reject the movement according to the
	  //// Metropolis-Boltzman criterion.
	  if (gsl_rng_uniform(gen) < exp(dE/T)){
		move(target,oldg,newg,part);
		part->module[target] = newg;
		E += dE;
	  }
	}// End of individual movements

	//// COLLECTIVE MOVEMENTS. ////
	for (i=collective_movements; i; collective_movements--){
	  //// MERGES //// 
	  //// Select two groups to merge
	  target = floor(gsl_rng_uniform(gen) * N);
	  g1 = part->module[target];
	  if (empty=-1){ // unless all nodes are together
		do{
		  target = floor(gsl_rng_uniform(gen) * (double)N);
		  g2 = part->module[target];
		} while (g1 == g2);


		// Compute the differences in energy.
		dE = 0.0;
		element = part->group[g1];
		while ((element = element->next) != NULL) {
		  element2 = part->group[g2];
		  while ((element2 = element2->next) != NULL) {
			dE += cmat[element->value][element2->value];
		  }
		}

		//// Accept or reject the movement according to the
		//// Metropolis-Boltzman criterion.
		if (gsl_rng_uniform(gen) < exp(dE/T)){
		  merge(g1,g2,part);
		  element = part->group[g2];
		  while ((element = element->next) != NULL)
			part->module[element->value] = g1;	  
		  E += dE;
		}
	   }// End unless all nodes are together (End of merge). 
	   //// SPLITS ////    
	   /* if (empty >= 0 ) { /\* if there are no empty groups, do nothing *\/ */
	   /* 	 /\* Select group to split *\/ */
	   /* 	 do { */
	   /* 	   target = floor(gsl_rng_uniform(gen) * (double)nnod); /\* node *\/ */
	   /* 	   target = nlist[target]->inGroup;    /\* target group *\/ */
	   /* 	 } while (part->group[target]->size == 1); */
		 
	   /* 	 // Split the group  */
	   /* 	 /\*GeneralSAGroupSplit(part->group[target], part->group[empty], */
	   /* 						 Ti, T, 0.95, */
	   /* 						 cluster_prob, */
	   /* 						 cmat, gen);*\/ */
		 
	   /* 	 // Calculate dE for splitting the groups. */
	   /* 	 dE = 0.0; */
	   /* 	 nod = part->group[target]; */
	   /* 	 while ((nod = nod->next) != NULL) { */
	   /* 	   nod2 = part->group[empty]; */
	   /* 	   while ((nod2 = nod2->next) != NULL) { */
	   /* 		 dE -= cmat[nod->node][nod2->node]; */
	   /* 	   } */
	   /* 	 } */
	   
	   /* 	 //// Accept or reject the movement according to the */
	   /* 	 //// Metropolis-Boltzman criterion.   */
	   /* 	 if ((dE > 0.0) && (gsl_rng_uniform(gen) > exp(dE/T))) */
	   /* 	   merge(dice,empty,part->group); // Revert the split. */
	   /* 	 else */
	   /* 	   E += dE; */
	   /* } // End of unless there is no empty group (=End of split). */
	}// End of collective movements

	//// BREAK THE LOOP IF NO CHANGES... ////
	if ( fabs(E - previousE) / fabs(previousE) < EPSILON_MOD
		 || fabs(previousE) < EPSILON_MOD)
	  {
	  nochange_count++;
	  // If we reach the limit...
	  if (nochange_count == nochange_limit){
		// If the current partition is the best so far. Terminate the
		// SA by breaking out of the temperature loop.
		if (E + EPSILON_MOD >= best_E) break;

		// Otherwise, reset the partition to the best one and proceed.
		FreeSimplePartition(part);
		part = CopySimplePartition(best_part);
		E = best_E;
		nochange_count = 0;
	  }
	  }
	// update the previous energy level.
	previousE = E;
	
	// Compare the current partition to the best partition so far and
	// update it if needed. 
	if ( E > best_E ) {
	  if (best_part!=NULL)
		FreeSimplePartition(best_part);
	  best_part = CopySimplePartition(part);

	  best_E = E;
	}

  } // End of the Temperature loop (end of SA). 

  grp = convertPartitionToGroup(part,net);
  FreeSimplePartition(best_part);
  FreeSimplePartition(part);
  return(grp);
}


  


/**
Given an array of linked lists, move the value val from the linked list old to new.
Return 0 if successfull. 
**/
unsigned int move(unsigned int val,
				  unsigned int old, unsigned int new,
				  struct partition *part){
  struct llist *element, *previous;
  
  element = part->group[old];
  previous = part->group[old];
  while (element->value != val && element->next != NULL){
	previous = element;
	element = element->next;
  }
  
  if (element->value != val)
	return(1);

  previous->next = element->next;
  element->next = part->group[new];
  part->group[new] = element;
  part->size[old]--;
  part->size[new]++;
  return(0);
}

/**
Merge two modules.
**/
void merge(unsigned int g1,
		   unsigned int g2,
		   struct partition *part){

  struct llist *element;
  unsigned int ishort, ilong;
  
  // Get the shorter list.
  if (part->size[g1]>part->size[g2]){
	ishort=g2;
	ilong=g1;
  }
  else{
	ishort=g1;
	ilong=g2;
  }
  
	
  // Traverse the shortest linked list.
  element = part->group[ishort];
  while ((element = element->next) != NULL){
	part->module[element->value] = ilong;
	part->size[ilong]++;
	part->size[ishort]--;
  }

  // Set the last element of the shorter list to the first of the longer list.
  element->next = part->group[ilong];

  // Set the shortest list at the place of the longest list (so the
  // values in part->module are correct)
  part->group[ilong] = part->group[ishort];
  
  // Remove the reference to the shortest list. 
  part->group[ishort] = NULL; 
}
