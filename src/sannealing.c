#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "bipartite.h"
#include "sannealing.h"
#include "costmatrix.h"

#define EPSILON_MOD 1.e-6

////////////////////////////////////////////////////////////////////////////////////:

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
	if (part->size[i] != 0){
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



struct partition * InitializePartition(unsigned int N){
  struct partition *part = NULL;
  unsigned int i;

  part = (struct partition *) malloc(sizeof(struct partition));
  part->nnod = N;
  part->nmod = N;
  part->nempty = 0;
  
  part->module = (unsigned int *) malloc(N*sizeof(unsigned int ));
  part->size = (unsigned int *) malloc(N*sizeof(unsigned int ));
  
  // Assign each node to a module.
  for(i=0;i<N;i++){
	part->module[i] = i;
	part->size[i] = 1;
  }
  return(part);
}

struct partition *CopySimplePartition(struct partition *part){
  struct partition *copy = NULL;
  unsigned int i, j;
  copy = (struct partition *) malloc(sizeof(struct partition));
  copy->nnod = part->nnod;
  copy->nmod = part->nmod;
  
  copy->module = (unsigned int *) malloc(part->nnod*sizeof(unsigned int ));
  copy->size = (unsigned int *) calloc(part->nnod,sizeof(unsigned int ));

  
  for(i=0;i<part->nnod;i++){
	// Copy the module number and update the size.
	j = part->module[i];
	copy->module[i] = j;
	copy->size[j] += 1;
  }
  return(copy);
}

void FreeSimplePartition(struct partition *part){
  free(part->module);
  free(part->size);
  free(part);
}
///////////////////////////////////////////////////////////////////////////////

/**
Given an array of linked lists, move the value val from the linked list old to new.
Return 0 if successfull. 
**/
unsigned int move(unsigned int target,
				  unsigned int old, unsigned int new,
				  struct partition *part){
  int s, sum, j;
  s = 0;
  for(j=0;j<part->nmod;j++){
	s+=part->size[j];
  }
  sum = s;
  if(sum != part->nnod){
	printf ("before: HHAAAAAAARG\n");
  }

  // Actually moving the node.
  part->module[target] = new;

  // Group size bookkeeping.
  if(!part->size[new]) part->nempty--;
  part->size[old]--;
  part->size[new]++;
  if(!part->size[old]) part->nempty++;

  s = 0;
  for(j=0;j<part->nmod;j++){
	s+=part->size[j];
  }
  sum = s;
  if(sum != part->nnod){
	printf ("after: HHAAAAAAARG\n");
  }

  return(0);
}

/**
Merge two modules.
**/
void merge(unsigned int g1,
		   unsigned int g2,
		   struct partition *part){

  unsigned int ishort, ilong, i;

  if (!part->size[g1]||!part->size[g2])
	return;
  
  // Get the shorter list.
  if (part->size[g1]>part->size[g2]){
	ishort=g2;
	ilong=g1;
  }
  else{
	ishort=g1;
	ilong=g2;
  }
  
  part->size[ilong] += part->size[ishort];
  for (i=0;i<part->nnod;i++){
	if (part->module[i] == ishort){
	  part->module[i] = ilong;
	  part->size[ishort]--;
	  if (part->size[ishort] == 0) break;
	}
  }
  part->nempty++;
  
}
double
dEAddNodeToGroup(unsigned int nodeid,
				 unsigned int groupid,
				 struct partition *part,
				 double **cmat){
  unsigned int i, count;
  double dE = 0.0;
  dE = 0.0;
  count = part->size[groupid];
  for (i = 0; i<part->nnod;i++){
	if(part->module[i]==groupid){
	  dE += cmat[nodeid][i];
	  count--;
	  if(count==0) break;
	}
  }
  return(dE);
}

double
dEMergeGroups(unsigned int group1,
			  unsigned int group2,
			  struct partition *part,
			  double **cmat){
  unsigned int i,j,count_long,count_short, gshort, glong;
  unsigned int *ishort, *ilong;
  double dE = 0.0;
  
  if (part->size[group1] > part->size[group2]){
	gshort = group2;
	glong = group1;
  }else{
	gshort = group1;
	glong = group2;
  }
  ishort = (unsigned int*) malloc(part->size[gshort]*sizeof(unsigned int));
  ilong = (unsigned int*) malloc(part->size[glong]*sizeof(unsigned int));

  count_short = part->size[gshort];
  count_long = part->size[glong];
  for (i = 0; i<part->nnod;i++){
	if(count_short && part->module[i]==gshort){
	  count_short--;
	  ishort[count_short] = i;
	  if(!count_short&&!count_long) break;}
	else if(count_long && part->module[i]==glong){
	  count_long--;
	  ilong[count_long] = i;
	  if(!count_short&&!count_long) break;}
	}
  
  dE = 0.0;
  for (i = 0; i<part->size[gshort];i++){
	for (j = 0; j<part->size[glong];j++){
	  dE += cmat[ishort[i]][ilong[j]];
	}
  }
  free(ishort);
  free(ilong);

  return(dE);
}



//////////////////////////////////////////////////////////////////////////////:

void
GeneralSAGroupSplit(unsigned int target, unsigned int empty,
					double Ti, double Tf, double Ts, double fac,
					double cluster_prob, unsigned int nochange_limit,
					double **cmat, struct partition *part,
					struct node_gra *net, gsl_rng *gen){
  unsigned int N, *indices, *modules, i, j, count, tnode;
  unsigned int oldg, newg, nochange_count = 0;
  double T, dE=0.0, E=0.0;
  N = part->size[target];
  indices = (unsigned int*) calloc(N,sizeof(unsigned int));
  modules = (unsigned int*) calloc(N,sizeof(unsigned int));

  
  count = N;
  j = 0;
  for (i=0; i<part->nnod; i++){
	if (part->module[i]==target){
	  indices[j];
	  count--;
	  j++;
	  if (!count) break;
	}
  }

  nochange_count = 0;
  for (T=Ti; T > Tf; T*=Ts) {
	//// INDIVIDUAL MOVEMENTS. ////
	  //// Select a node.
	  tnode = floor(gsl_rng_uniform(gen) * (double)N);
	  tnode = indices[tnode];

	  // Switch its group.
	  oldg = part->module[tnode];
	  if (oldg == target)
		newg = empty;
	  else
		newg = target;
	  
	  //// Computing the difference in energy.
	  // Removing the node from its group.
	  dE = -dEAddNodeToGroup(tnode,oldg,part,cmat);
	  // Adding it to its new group.
	  dE += dEAddNodeToGroup(tnode,newg,part,cmat);

	  //// Accept or reject the movement according to the
	  //// Metropolis-Boltzman criterion.
	  if (gsl_rng_uniform(gen) < exp(dE/T)){
		move(tnode,oldg,newg,part);
		E += dE;
	  } else{
		dE =0;
	  }
	  if (fabs(dE) / fabs(E) < EPSILON_MOD
		  || fabs(E) < EPSILON_MOD){
		nochange_count++;
		if(nochange_count>nochange_limit){
		  break;
		}
	  }
  }// End of SA.
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
  unsigned int target, oldg, newg, i, g1, g2, j;
  int empty=0;
  double T = Ti;
  double cluster_prob = 0.5;
  unsigned int ngroup = N, nochange_count=0, nochange_limit = 25;
  double dE=0.0, E=0.0, previousE=0.0, best_E=-1.0/0.0; //initial best is -infinity.
  struct llist *element, *element2;
  part = InitializePartition(N);
  int s=0;
  
  // Get the number of individual and collective movements. 
  iterationNumber(N, fac, &individual_movements, &collective_movements);
  
  /// SIMULATED ANNEALING ///
  //printf ("#T\tE\n");
  for (T=Ti; T > Tf; T = T*Ts) {
	//printf ("%1.5e\t%1.5e \n",T,E);
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
	  dE = -dEAddNodeToGroup(target,oldg,part,cmat);
	  // Adding it to its new group.
	  dE += dEAddNodeToGroup(target,newg,part,cmat);

	  //// Accept or reject the movement according to the
	  //// Metropolis-Boltzman criterion.
	  if (gsl_rng_uniform(gen) < exp(dE/T)){
		move(target,oldg,newg,part);
		E += dE;
	  }
	}// End of individual movements

	//// COLLECTIVE MOVEMENTS. ////
	for (i=collective_movements; i; i--){
	  //// MERGES //// 
	  //// Select two groups to merge
	  target = floor(gsl_rng_uniform(gen) * N);
	  g1 = part->module[target];
	  if (part->nempty<part->nmod-1){ // unless all nodes are together
		do{
		  target = floor(gsl_rng_uniform(gen) * (double)N);
		  g2 = part->module[target];
		} while (g1 == g2);


		// Compute the differences in energy.
		dE = dEMergeGroups(g1,g2,part,cmat);

		//// Accept or reject the movement according to the
		//// Metropolis-Boltzman criterion.
		if (gsl_rng_uniform(gen) < exp(dE/T)){
		  merge(g1,g2,part);
		  E += dE;
		}
	   }// End unless all nodes are together (End of merge). 
	   //// SPLITS ////    
	  if (part->nempty) { //unless there is no empty groups.
		// Select the empty group;
		
		for (j=0;j<part->nmod;j++){
		  if(part->size[j]==0)
			empty = j;
		}
		
	   	 do {
	   	   target = floor(gsl_rng_uniform(gen) * (double)part->nnod); //Node
	   	   target = part->module[target]; // Group. 
	   	 } while (part->size[target] == 1);
		 
	   	 // Split the group
	   	 GeneralSAGroupSplit(target, empty,
	   						 Ti, Tf, 0.95, fac,
	   						 cluster_prob, nochange_limit,
	   						 cmat, part, net, gen);
		 
	   	 // Calculate dE for splitting the groups.
	   	 dE = -dEMergeGroups(target,empty,part,cmat);
	   
	   	 //// Accept or reject the movement according to the
	   	 //// Metropolis-Boltzman criterion.
	   	 if ((dE > 0.0) && (gsl_rng_uniform(gen) < exp(dE/T)))
	   	   E += dE;
	   	 else
		   merge(target,empty,part); // Revert the split.
	   } // End of unless there is no empty group (=End of split).
	} //End of collective movements

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
  printf ("Partition\n");
  for (i=0;i<part->nnod;i++){
	printf (" %d ",part->module[i]);
  }
  printf ("\n");
  grp = convertPartitionToGroup(part,net);
  FreeSimplePartition(best_part);
  FreeSimplePartition(part);
  return(grp);
}

//////////////////////////////////////////////////////::

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
  int i, j;
  for(i=0;i<N;i++){
	for(j=0;j<N;j++){
	  printf ("\t%e\t",cmat[i][j]);	  
	}
	printf ("\n");
  }
  
  ret = GeneralSA(cmat, N, net,
  				  fac, Ti, Tf,  Ts,
  				  gen);
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
  int i,j;
   for(i=0;i<N;i++){
	for(j=0;j<N;j++){
	  printf ("\t%e\t",cmat[i][j]);	  
	}
	printf ("\n");
  }
 
  ret = GeneralSA(cmat, N, net,
  				  fac, Ti, Tf,  Ts,
  				  gen);
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
  int i,j;
   for(i=0;i<N;i++){
	for(j=0;j<N;j++){
	  printf ("\t%e\t",cmat[i][j]);	  
	}
	printf ("\n");
  }
 
  projected = ProjectBipart(binet);
  ret = GeneralSA(cmat, N, projected,
  				  fac, Ti, Tf,  Ts,
  				  gen);
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
  int i,j;
   for(i=0;i<N;i++){
	for(j=0;j<N;j++){
	  printf ("\t%e\t",cmat[i][j]);	  
	}
	printf ("\n");
  }
 
  projected = ProjectBipart(binet);
  ret = GeneralSA(cmat, N, projected,
  				  fac, Ti, Tf,  Ts,
  				  gen);
  free(cmat);
  return(ret);
}
