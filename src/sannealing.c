#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#include "sannealing.h"
#include "partition.h"
#include "movements.h"

#define VERBOSE
#ifndef EXPLAIN
#define explain(M, ...)
#else
#define explain(M, ...) fprintf(stderr,M, ##__VA_ARGS__)
#endif

#ifndef VERBOSE
#define info(M, ...)
#else
#define info(M, ...) fprintf(stderr,M, ##__VA_ARGS__)
#endif

#define EPSILON_MOD 1.e-6

/**
Modularity optimisation using simulated annealing.

@param part Initial partition.
@param adj Adjacency array.
@param fac Iteration factor.
@param Ti Initial temperature.
@param Tf Finial temperature.
@param Ts Cooling factor.
@param proba_components probability to try using connected components to split a module.
@param nochange_limit number of consecutive non improving step before stopping.
@param gen random number generator.
**/
unsigned int
GeneralSA(Partition **ppart, AdjaArray *adj,
		  double fac,
		  double Ti, double Tf, double Ts,
		  double proba_components,
		  unsigned int nochange_limit,
		  gsl_rng *gen)
{
  Partition *part = *ppart;
  Partition *best_part = NULL;
  unsigned int individual_movements, collective_movements;
  unsigned int target, oldg, newg, i, g1, g2, j, empty, ncomponent, err;
  double T = Ti;
  unsigned int nochange_count=0;

  double dE=0.0, E=0.0, previousE=0.0, best_E=-1.0/0.0;  //initial best is -infinity.
  info ("#Simulated annealing:\n");
  explain ("#Ti: %f Tf: %f c: %f fac: %f\n",Ti,Tf,Ts,fac);
  explain ("#nochang_limit: %d, proba_components: %f \n",nochange_limit, proba_components);

  // Get the number of individual and collective movements.
  if (fac * (double)(part->N * part->N) < 10)
	individual_movements = 10;
  else
	individual_movements = floor(fac * (double)(part->N * part->N));
  if (fac * (double)part->N < 2)
	collective_movements = 2;
  else
	collective_movements = floor(fac * (double)part->N);

  info ("#T\tE\tStop\n");
  /// SIMULATED ANNEALING ///
  for (T=Ti; T > Tf; T = T*Ts) {
	info ("%e\t%e\t%d\n", T , E, nochange_count);
	//// INDIVIDUAL MOVEMENTS. ////
	for (i=individual_movements; i; i--) {
	  //// Select a node and a target group.
	  target = floor(gsl_rng_uniform(gen) * (double)part->N);
	  oldg = part->nodes[target]->module;
	  do {
		newg = floor(gsl_rng_uniform(gen) * part->M);
	  } while (newg == oldg);

	  //// Computing the difference in energy.
	  dE = dEChangeModule(target,newg,part,adj);

	  //explain("Moving node %d from %d to %d, %e\n",target,oldg,newg,dE);
	  //// Accept or reject the movement according to the
	  //// Metropolis-Boltzman criterion.
	  if ( (dE>=0) || (gsl_rng_uniform(gen) < exp(dE/T)) ){
		ChangeModule(target,newg,part);
		E += dE;
	  }
	}// End of individual movements

	//// COLLECTIVE MOVEMENTS. ////

	for (i=collective_movements; i; i--){
	  //// MERGES ////
	  //// Select two groups to merge
	  target = floor(gsl_rng_uniform(gen) * part->N);
	  g1 = part->nodes[target]->module;
	  if (part->nempty<part->M-1){ // unless all nodes are together
		do{
		  target = floor(gsl_rng_uniform(gen) * (double)part->N);
		  g2 = part->nodes[target]->module;
		} while (g1 == g2);


		// Compute the differences in energy.
		dE = dEMergeModules(g1,g2,part,adj);
		explain("Merging module %d (%d) and %d (%d), %e\n",g1,part->modules[g1]->size,g2,part->modules[g2]->size,dE);
		//// Accept or reject the movement according to the
		//// Metropolis-Boltzman criterion.
		if ( (dE>=0) || (gsl_rng_uniform(gen) < exp(dE/T)) ){
		  MergeModules(g1,g2,part);
		  E += dE;
		}
	  }// End unless all nodes are together (End of merge).

	   //// SPLITS ////
	  if (part->nempty) { //unless there is no empty groups.
		// Select the first empty group;
		for (j=0;j<part->M;j++){
		  if(!part->modules[j]->size){
			empty = j;
			break;
		  }
		}

		// select the group to split.
		do {
		  target = floor(gsl_rng_uniform(gen) * (double)part->N); //Node
		  target = part->nodes[target]->module; // Group.
		} while (part->modules[target]->size == 1);

		// Split the group

		// Try to split the module along connect components with a
		// given probability.
		explain("Split modules %d (%d) and %d (%d) Proba: (%f)\n",target, part->modules[target]->size, empty,part->modules[empty]->size,proba_components);
		ncomponent = 1;
		if (gsl_rng_uniform(gen)<proba_components){
		  ncomponent = SplitModuleByComponent(target, empty,
											  part,adj,
											  gen);
		  explain("Using connected components: %d cc\n",ncomponent);
		}

		// If the previous split was unsuccessful (or not even tried)
		// try to split it with a nested simulated annealing. Start
		// with a high temperature (Ti) and do individual movements
		// until it reaches the global SA temperature (T).
		if (ncomponent==1){
		  err = SplitModuleSA(target, empty,
							  Ti, T, 0.95,
							  nochange_limit,
							  part, adj, gen);
		  if (err)
			return 2;
		  explain("Using a nested SA.\n");
		}

		// If the split actually gave two non empty modules...
		if ((part->modules[target]->size) && (part->modules[empty]->size)){

		  // The energy variation when splitting a group is equal to the
		  // opposite of the energy variation when merging them back.
		  dE = -dEMergeModules(target,empty,part,adj);

		  explain("Splitted modules: %d (%d) and %d (%d), %e",
				  target, part->modules[target]->size,
				  empty, part->modules[empty]->size,
				  dE);

		  // Accept or revert the movement according to the
		  // Metropolis-Boltzman criterion.
		  if ((dE >= 0) || (gsl_rng_uniform(gen) < exp(dE/T))){
			E += dE;
			explain("... Accepted split\n");
		  } else{
			MergeModules(target,empty,part); // Revert the split.
			explain("... Reverted split\n");
		  }
		} else { //If one of the module post split is still empty.
		  explain("Trivial split\n");
		}

	  } // End of unless there is no empty group (=End of split).
	} //End of collective movements

	//// BREAK THE LOOP IF NO CHANGES... ////
	if ( (fabs(E - previousE) / fabs(previousE) < EPSILON_MOD)
		 || (fabs(previousE) < EPSILON_MOD) )
	  {
		nochange_count++;
		// If we reach the limit...
		if (nochange_count == nochange_limit){
		  info ("# Too much rounds without changes (%d)... \n",nochange_count);
		  // If the current partition is the best so far. Terminate the
		  // SA by breaking out of the temperature loop.
		  if (E+EPSILON_MOD>= best_E) break;

		  // Otherwise, reset the partition to the best one and proceed.
		  info ("# Restarting from a better place (%e<%e)\n",E,best_E);
		  E = best_E;
		  nochange_count = 0;
		  FreePartition(part);
		  part = CopyPartitionStruct(best_part);
		}
	  }
	// Compare the current partition to the best partition so far and
	// update it if needed.
	else if ( E > best_E) {
	  info ("# Saving a new best partition (%e)\n",E);
	  if (best_part!=NULL)
		FreePartition(best_part);
	  best_part = CopyPartitionStruct(part);
	  best_E = E;
	  nochange_count = 0;
	}
   	// update the previous energy level.
	previousE = E;

  } // End of the Temperature loop (end of SA).

  info ("# End of SA, best partition so far: %e\n",best_E);

  FreePartition(part);
  *ppart = best_part;
  return(0);
}



/**
Split the target module into two using a small nested simulated
annealing.

This function start by randomly spliting the target modules into two
modules (using the empty module id provided) and then perform an
individual movement-only simulated annealing to get the best possible
split.

Note that the modularity cost are computed using the values in part
and adj, so they are related to the global SA.
**/
unsigned int
SplitModuleSA(unsigned int target, unsigned int empty,
			  double Ti, double Tf, double Ts,
			  unsigned int nochange_limit,
			  Partition *part, AdjaArray *adj,
			  gsl_rng *gen){

  unsigned int N, *indices, nodeid, i;
  unsigned int oldg, newg, nochange_count = 0;
  Node * node;
  double T, dE=0.0;

  N = part->modules[target]->size;

  // Build an array of indices to be able to draw one node at random.
  indices = (unsigned int*) calloc(N,sizeof(unsigned int));
  if (indices == NULL){
	perror("Error while splitting module");
	return 1;}
  for(node=part->modules[target]->first,  i = 0; node!=NULL; node = node->next,i++)
	indices[i] = node->id;

  // Split the module randomly into two.
  for(i = 0; i<N; i++){
	if (gsl_rng_uniform(gen) < 0.5)
	   ChangeModule(indices[i],empty,part);
  }

  for (T=Ti; T > Tf; T*=Ts) {
	//// Select a random node in the module.
	nodeid = floor(gsl_rng_uniform(gen) * (double)N);
	nodeid = indices[nodeid];

	// Switch its group.
	oldg = part->nodes[nodeid]->module;
	newg = oldg == target? empty : target;

	//// Computing the difference in energy.
	dE = dEChangeModule(nodeid,newg,part,adj);

	//// Accept or reject the movement according to the
	//// Metropolis-Boltzman criterion.
	if ((dE>=0) || (gsl_rng_uniform(gen) < exp(dE/T)))
	  ChangeModule(nodeid,newg,part);
	else
	  dE = 0;

	// If the change was to small, update the nochange count and break
	// out the loop if the limit is reached.
	if (fabs(dE) < EPSILON_MOD){
	  nochange_count++;
	  if (nochange_count>nochange_limit) break;
	} else{
	  nochange_count=0;
	}
  }// End of SA.
  free(indices);
  return 0;
}
