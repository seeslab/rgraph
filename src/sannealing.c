#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "bipartite.h"

#define EPSILON_MOD 1.e-6
#define EPSILON_MOD_B 1.e-6


/*
  ---------------------------------------------------------------------
  Given a group, SAGroupSplit returns a split of the nodes into two
  subgroups. It uses SA, so initial and final temperature must be
  provided. If cluster_sw == 1, then the function checks first whether
  the group is disconnected; if it is, then it returns (with some
  probability) a partition along the lines of the existing clusters
  without any SA (with complementary probability, it still does the
  SA).
  ---------------------------------------------------------------------
*/
struct group *
SAGroupSplit(struct group *targ,
			   double Ti, double Tf, double Ts,
			   int cluster_sw,
			   gsl_rng *gen)
{
  struct group *glist[2];
  struct group *split = NULL;
  struct node_gra **nodeList;
  struct node_gra *net = NULL;
  struct node_gra *p = NULL;
  struct node_p;
  int nnod = 0;
  int i;
  int des;
  int target, oldg, newg;
  int innew, inold, nlink;
  int totallinks=0;
  double dE=0.0, energy=0.0;
  double T;
  int ngroups, g1, g2;
  double prob = 0.5;

  nodeList = (struct node_gra**)calloc(targ->size,
					   sizeof(struct node_gra *));
  glist[0] = NULL;
  glist[1] = NULL;

  /* Build a network from the nodes in the target group */
  net = BuildNetFromGroup(targ);

  /* Check if the network is connected */
  split = ClustersPartition(net);
  ngroups = NGroups(split);

  if (
	  cluster_sw == 1 &&         /*  cluster switch is on */
	  ngroups > 1 &&             /*  network is not connected */
	  gsl_rng_uniform(gen) < prob  /*  with some probability */
	  ) {

	/* Merge groups randomly until only two are left */
	while (ngroups > 2) {
	  /* Select two random groups */
	  g1 = ceil(gsl_rng_uniform(gen)* (double)ngroups);
	  do {
	g2 = ceil(gsl_rng_uniform(gen)* (double)ngroups);
	  } while (g2 == g1);

	  glist[0] = split;
	  for(i=0; i<g1; i++)
	glist[0] = glist[0]->next;
	  glist[1] = split;
	  for(i=0; i<g2; i++)
	glist[1] = glist[1]->next;

	  /* Merge */
	  MergeGroups(glist[0], glist[1]);
	  split = CompressPart(split);
	  ngroups--;
	}
  }

  else { /*  Network is connected or we want to ignore the clusters */
	/* Remove SCS partition */
	RemovePartition(split);
	ResetNetGroup(net);

	/* Create the groups */
	split = CreateHeaderGroup();
	glist[0] = CreateGroup(split, 0);
	glist[1] = CreateGroup(split, 1);

	/* Randomly assign the nodes to the groups */
	p = net;
	while ((p = p->next) != NULL) {
	  nodeList[nnod] = p;
	  totallinks += CountLinks(p);
	  nnod++;

	  des = floor(gsl_rng_uniform(gen) * 2.0);
	  AddNodeToGroup(glist[des], p);
	}
	totallinks /= 2;

	/* Do the SA to "optimize" the splitting */
	if (totallinks > 0) {
	  T = Ti;
	  while (T > Tf) {

	for (i=0; i<nnod; i++) {
	  target = floor(gsl_rng_uniform(gen) * (double)nnod);
	  oldg = nodeList[target]->inGroup;
	  if(oldg == 0)
		newg = 1;
	  else
		newg = 0;

	  /* Calculate the change of energy */
	  inold = NLinksToGroup(nodeList[target],glist[oldg]);
	  innew = NLinksToGroup(nodeList[target],glist[newg]);
	  nlink = CountLinks(nodeList[target]);

	  dE = 0.0;
	  dE -= (double)(2 * glist[oldg]->inlinks) /
		(double)totallinks -
		(double)(glist[oldg]->totlinks+glist[oldg]->inlinks) *
		(double)(glist[oldg]->totlinks+glist[oldg]->inlinks) /
		((double)totallinks * (double)totallinks);
	  dE -= (double)(2 * glist[newg]->inlinks) /
	  (double)totallinks -
		(double)(glist[newg]->totlinks+glist[newg]->inlinks) *
		(double)(glist[newg]->totlinks+glist[newg]->inlinks) /
		((double)totallinks * (double)totallinks);
	  dE += (double)(2*glist[oldg]->inlinks - 2*inold) /
		(double)totallinks -
		(double)(glist[oldg]->totlinks + glist[oldg]->inlinks -
			 nlink ) *
		(double)(glist[oldg]->totlinks + glist[oldg]->inlinks -
			 nlink ) /
		((double)totallinks * (double)totallinks);
	  dE += (double)(2*glist[newg]->inlinks + 2*innew) /
		(double)totallinks -
		(double)(glist[newg]->totlinks + glist[newg]->inlinks +
			 nlink ) *
		(double)(glist[newg]->totlinks + glist[newg]->inlinks +
			 nlink ) /
		((double)totallinks * (double)totallinks);

	  /* Accept the move according to Metropolis */
	  if( (dE >=0.0) || (gsl_rng_uniform(gen) < exp(dE/T)) ){
		MoveNode(nodeList[target],glist[oldg],glist[newg]);
		energy += dE;
	  }
	}

	T = T * Ts;
	  } /*  End of temperature loop */
	} /*  End if totallinks > 0 */
  }

  /* Free memory */
  RemoveGraph(net);
  free(nodeList);

  /* Done */
  return split;
}


/*
  ---------------------------------------------------------------------
  Given a network, use simulated annealing to find the partition that
  maximizes the Newman-Girvan modularity. Ti, Tf, and Ts are the
  initial and final temperatures, and Ts is the temperature cooling
  factor. ngroup is the number of modules used (0=use the number of
  nodes in the network). initial_sw defines the initial configuration
  of nodes into modules ('o'=one node per module, 'r'=random placemnt
  of nodes into modules). If initial_sw is set to 'o', then ngroup is
  overridden and set to the number of nodes in the
  network. collective_sw determines whether collective merge-split
  moves are used (1) or not (0). output_sw determines the amount of
  information printed to stdout (from less information to more
  information: 'n'=none, 'b'=backup, 'm'=minimal, 's'=backup +
  minimal, 'v'=verbose, and 'd'=debug).
  ---------------------------------------------------------------------
*/
struct group *
SACommunityIdent(struct node_gra *net,
		 double Ti, double Tf, double Ts,
		 double fac,
		 int ngroup,
		 char initial_sw,
		 int collective_sw,
		 char output_sw,
		 gsl_rng *gen)
{
  int i;
  struct group *part = NULL;
  struct group *split = NULL, *g = NULL;
  struct group **glist = NULL, *lastg;
  struct node_gra **nlist;
  struct node_gra *p;
  struct node_lis *nod;
  int dice;
  int empty;
  int newg, oldg;
  int nnod;
  int totallinks = 0;
  int innew,inold,nlink;
  double energy = -1.0, dE = 0.0;
  double T;
  int g1, g2;
  double energyant = 0.0;
  int count = 0, limit = 25; /*  to stop the search if the energy does
				 not change */
  int cicle1, cicle2;
  struct group *best_part = NULL;
  double best_E = -100.0;
  void *nodeDict;
  FILE *outf;

  /*
	Preliminaries: Initialize, allocate memory, and place nodes in
	initial groups
	-------------------------------------------------------------------
  */
  /* Initialize some variables */
  nnod = CountNodes(net);
  ResetNetGroup(net);
  part = CreateHeaderGroup();
  p = net;

  /* Create a node dictionary for fast access to nodes by label */
  nodeDict = MakeLabelDict(net);

  /* Allocate memory for the node list */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));

  /* Create the groups and assign nodes to the initial group according
	 to initial_sw. Additionally, map nodes and groups to lists for
	 faster access, and count the total number of links. */
  switch (initial_sw) {

  case 'o':         /*  One node in each group */
	ngroup = nnod;
	glist = (struct group **) calloc(ngroup, sizeof(struct group *));
	while ((p = p->next) != NULL) {
	  glist[p->num] = CreateGroup(part, p->num);
	  nlist[p->num] = p;
	  AddNodeToGroup(glist[p->num], p);
	  totallinks += CountLinks(p);
	}
	break;

  case 'r':        /*  Random placement of nodes in groups */
	glist = (struct group **) calloc(ngroup, sizeof(struct group *));
	lastg = part;
	for (i=0; i<ngroup; i++)
	  glist[i] = lastg = CreateGroup(lastg, i);
	while ((p = p->next) != NULL) {
	  nlist[p->num] = p;
	  dice = floor(gsl_rng_uniform(gen)* (double)ngroup);
	  AddNodeToGroup(glist[dice], p);
	  totallinks += CountLinks(p);
	}
	break;
  }

  /*
	Determine the number of iterations at each temperature
	-------------------------------------------------------------------
  */
  if (fac * (double)(nnod * nnod) < 10)
	cicle1 = 10;
  else
	cicle1 = floor(fac * (double)(nnod * nnod));

  if (fac * (double)nnod < 2)
	cicle2 = 2;
  else
	cicle2 = floor(fac * (double)nnod);

  /*
	Do the simulated annealing
	-------------------------------------------------------------------
  */
  /* Determine initial values */
  T = Ti;
  energy = Modularity(part);

  /* Temperature loop */
  while (T > Tf && count < limit) {

	/* Output */
	switch (output_sw) {
	case 'n':
	  break;
	case 'b':
	  break;
	case 's':
	  fprintf(stderr, "%g %lf %g\n",1.0/T, energy, T);
	  break;
	case 'm':
	  fprintf(stderr, "%g %lf %g\n",1.0/T, energy, T);
	  break;
	case 'v':
	  fprintf(stderr, "%g %lf %lf %g\n",
		  1.0/T, energy, Modularity(part), T);
	  break;
	case 'd':
	  FPrintPartition(stderr, part, 0);
	  fprintf(stderr, "%g %lf %lf %g\n",
		  1.0/T, energy, Modularity(part), T);
	}

	/* Do all the individual moves */
	for (i=0; i<cicle1; i++) {

	  /* Propose an individual move */
	  dice = floor(gsl_rng_uniform(gen) * (double)nnod);
	  oldg = nlist[dice]->inGroup;
	  do {
	newg = floor(gsl_rng_uniform(gen) * (double)ngroup);
	  } while (newg == oldg);

	  /* Calculate the change of energy */
	  inold = NLinksToGroup(nlist[dice], glist[oldg]);
	  innew = NLinksToGroup(nlist[dice], glist[newg]);
	  nlink = CountLinks(nlist[dice]);
	  dE = 0.0;
	  dE -= (double)(2 * glist[oldg]->inlinks) /
	(double)totallinks -
	(double)(glist[oldg]->totlinks+glist[oldg]->inlinks) *
	(double)(glist[oldg]->totlinks+glist[oldg]->inlinks) /
	((double)totallinks * (double)totallinks);
	  dE -= (double)(2 * glist[newg]->inlinks) /
	(double)totallinks -
	(double)(glist[newg]->totlinks+glist[newg]->inlinks) *
	(double)(glist[newg]->totlinks+glist[newg]->inlinks) /
	((double)totallinks * (double)totallinks);
	  dE += (double)(2*glist[oldg]->inlinks - 2*inold) /
	(double)totallinks -
	(double)(glist[oldg]->totlinks + glist[oldg]->inlinks -
		 nlink ) *
	(double)(glist[oldg]->totlinks + glist[oldg]->inlinks -
		 nlink ) /
	((double)totallinks * (double)totallinks);
	  dE += (double)(2*glist[newg]->inlinks + 2*innew) /
	(double)totallinks -
	(double)(glist[newg]->totlinks + glist[newg]->inlinks +
		 nlink ) *
	(double)(glist[newg]->totlinks + glist[newg]->inlinks +
		 nlink ) /
	((double)totallinks * (double)totallinks);

	  /* Accept the change according to Metroppolis */
	  if ((dE >= 0.0) || (gsl_rng_uniform(gen) < exp(dE/T))) {
	MoveNode(nlist[dice], glist[oldg], glist[newg]);
	energy += dE;
	  }
	} /*  End of individual moves */

	/* Do all the collective moves */
	if (collective_sw == 1) {
	  for ( i=0; i < cicle2; i++ ){

	/* MERGE */
	dice = floor(gsl_rng_uniform(gen) * (double)nnod);
	g1 = nlist[dice]->inGroup;

	if (glist[g1]->size < nnod) { /*  Unless all nodes are together */
	  do{
		dice = floor(gsl_rng_uniform(gen) * (double)nnod);
		g2 = nlist[dice]->inGroup;
	  } while (g1 == g2);

	  /* Calculate the change of energy */
	  nlink = NG2GLinks(glist[g1], glist[g2]);
	  dE = 0.0;
	  dE -= (double)(2*glist[g1]->inlinks) / (double)totallinks -
		((double)(glist[g1]->totlinks + glist[g1]->inlinks) *
		 (double)(glist[g1]->totlinks + glist[g1]->inlinks) ) /
		((double)totallinks * (double)totallinks );
	  dE -= (double)(2*glist[g2]->inlinks) / (double)totallinks -
		((double)(glist[g2]->totlinks + glist[g2]->inlinks) *
		 (double)(glist[g2]->totlinks + glist[g2]->inlinks) ) /
		((double)totallinks * (double)totallinks );
	  dE += 2.0*(double)(glist[g1]->inlinks +
				 glist[g2]->inlinks+nlink)
		/ (double)totallinks -
		(double)(glist[g1]->totlinks + glist[g1]->inlinks +
			 glist[g2]->totlinks + glist[g2]->inlinks ) *
		(double)(glist[g1]->totlinks + glist[g1]->inlinks +
			 glist[g2]->totlinks + glist[g2]->inlinks ) /
		((double)totallinks * (double)totallinks);

	  /* Accept the change according to Metroppolis */
	  if ((dE >= 0.0) || (gsl_rng_uniform(gen) < exp(dE/T))) {
		MergeGroups(glist[g1], glist[g2]);
		energy += dE;
	  }
	}

	/* SPLIT */
	dice = floor(gsl_rng_uniform(gen) * (double)nnod); /*  target node */
	dice = nlist[dice]->inGroup;                     /*  target group */

	/* Look for an empty group */
	g = part;
	empty = -1;
	while (((g = g->next) != NULL) && (empty < 0))
	  if (g->size == 0)
		empty = g->label;

	if (empty >= 0) { /*  if there are no empty groups, do nothing */
	  /* Find a reasonable split */
	  split = SAGroupSplit(glist[dice], Ti, T, 0.95, 1, gen);

	  /* Split the group */
	  nod = (split->next)->nodeList;
	  while ((nod = nod->next) != NULL) {
		MoveNode(GetNodeDict(nod->nodeLabel, nodeDict),
			 glist[dice],
			 glist[empty]);
	  }
	  RemovePartition(split);
	  split = NULL;

	  /* Calculate the change of energy associated to remerging
		 the groups */
	  nlink = NG2GLinks(glist[dice], glist[empty]);
	  dE = 0.0;
	  dE -= (double)(2*glist[dice]->inlinks) /
		(double)totallinks -
		((double)(glist[dice]->totlinks+glist[dice]->inlinks) *
		 (double)(glist[dice]->totlinks+glist[dice]->inlinks))/
		((double)totallinks * (double)totallinks );
	  dE -= (double)(2*glist[empty]->inlinks) /
		(double)totallinks -
		((double)(glist[empty]->totlinks + glist[empty]->inlinks) *
		 (double)(glist[empty]->totlinks + glist[empty]->inlinks))/
		((double)totallinks * (double)totallinks );
	  dE += 2.0*(double)(glist[dice]->inlinks+
				 glist[empty]->inlinks+nlink)
		/ (double)totallinks -
		(double)(glist[dice]->totlinks + glist[dice]->inlinks
			 + glist[empty]->totlinks +
			 glist[empty]->inlinks ) *
		(double)(glist[dice]->totlinks + glist[dice]->inlinks +
			 glist[empty]->totlinks + glist[empty]->inlinks ) /
		((double)totallinks * (double)totallinks);

	  /* Accept the change according to "inverse" Metroppolis.
		 Inverse means that the algorithm is applied to the split
		 and NOT to the merge! */
	  if ((dE > 0.0) && (gsl_rng_uniform(gen) > exp(-dE/T))) {
		MergeGroups(glist[dice], glist[empty]);
	  }
	  else{
		energy -= dE;
	  }
	}  /*  End of if empty >= 0 */
	  }  /*  End of collective moves */
	}  /*  End of if collective_sw==1 */

	/* Update the no-change counter */
	if (fabs(energy - energyant) / fabs(energyant) < EPSILON_MOD ||
	fabs(energyant) < EPSILON_MOD) {
	  count++;

	  /* If the SA is ready to stop (count==limit) but the current
	 partition is not the best one so far, replace the current
	 partition by the best one and continue from there. */
	  if ((count == limit) && (energy + EPSILON_MOD < best_E)) {
	switch (output_sw) {
	case 'n':
	  break;
	case 'b':
	  break;
	default:
	  fprintf(stderr, "# Resetting partition\n");
	  break;
	}

	/* Remap the partition to the best partition */
	RemovePartition(part);
	g = part = CopyPartition(best_part);
	while ((g = g->next) != NULL)
	  glist[g->label] = g;
	MapPartToNet(part, net);

	/* Reset energy and counter */
	energy = best_E;
	count = 0;
	  }
	}

	else {
	  count = 0;
	}
	/* Update the last energy */
	energyant = energy;

	/* Compare the current partition to the best partition so far and
	   save the current if it is better than the best so far. */
	if (energy > best_E) {
	  if (best_part != NULL)
	RemovePartition(best_part);
	  best_part = CopyPartition(part);
	  MapPartToNet(part, net); /*  MUST DO this after copying a
				   part! */
	  best_E = energy;
	}

	/* Save the partition to a file if necessary */
	switch (output_sw) {
	case 'b':
	  outf = fopen("part.tmp", "w");
	  FPrintPartition(outf, best_part, 1);
	  fclose(outf);
	case 's':
	  outf = fopen("part.tmp", "w");
	  FPrintPartition(outf, best_part, 1);
	  fclose(outf);
	default:
	  break;
	}

	/* Uptade the temperature */
	T = T * Ts;

  } /*  End of simulated annealing */

  /* Free memory */
  RemovePartition(best_part);
  FreeLabelDict(nodeDict);
  free(glist);
  free(nlist);

  /* Done */
  return CompressPart(part);
}


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"

#define EPSILON_MOD 1.e-6

struct group *ThermalPercNetworkSplitWeight(struct group *targ,
						double Ti, double Tf,
						gsl_rng *gen)
{
  struct group *glist[2];
  struct group *split = NULL;
  struct node_gra **nlist;
  struct node_gra *net = NULL;
  struct node_gra *p = NULL;
  struct node_p;
  int nnod = 0;
  int i;
  int des;
  int target,oldg,newg;
  double innew,inold,nlink;
  double totallinks=0.0;
  double dE=0.0, energy=0.0;
  double T, Ts = 0.95;
  int ngroups, g1, g2;
  double prob = 0.5;

  nlist = (struct node_gra**)calloc(targ->size,
					   sizeof(struct node_gra *));

  glist[0] = NULL;
  glist[1] = NULL;

  // Build a network from the nodes in the target group
  net = BuildNetFromGroup(targ);

  // Check if the network is connected
  split = ClustersPartition(net);
  ngroups = NGroups(split);

  if ( ngroups > 1 && gsl_rng_uniform(gen) < prob) { // Network is not
						   // connected

	// Merge groups randomly until only two are left
	while (ngroups > 2) {
	  // Select two random groups
	  g1 = ceil(gsl_rng_uniform(gen)* (double)ngroups);
	  do {
	g2 = ceil(gsl_rng_uniform(gen)* (double)ngroups);
	  } while (g2 == g1);

	  glist[0] = split;
	  for(i=0; i<g1; i++)
	glist[0] = glist[0]->next;
	  glist[1] = split;
	  for(i=0; i<g2; i++)
	glist[1] = glist[1]->next;

	  // Merge
	  MergeGroups(glist[0], glist[1]);
	  split = CompressPart(split);
	  ngroups--;
	}
  }

  else { // Network IS connected
	// Remove SCS partition
	RemovePartition(split);
	ResetNetGroup(net);

	// Create the groups
	split = CreateHeaderGroup();
	glist[0] = CreateGroup(split,0);
	glist[1] = CreateGroup(split,1);

	// Randomly assign the nodes to the groups
	p = net;
	while(p->next != NULL){
	  p = p->next;
	  nlist[nnod] = p;
	  totallinks += NodeStrength(p);
	  nnod++;

	  des = floor(gsl_rng_uniform(gen)*2.0);
	  AddNodeToGroup(glist[des],p);
	}

	totallinks /= 2.0;

	// Do the SA to "optimize" the splitting
	if ( totallinks > 0 ) {
	  T = Ti;
	  while( T > Tf){

	for (i=0; i< nnod; i++){
	  target = floor(gsl_rng_uniform(gen) * (double)nnod);
	  oldg = nlist[target]->inGroup;
	  if(oldg == 0)
		newg = 1;
	  else
		newg = 0;

	  // Calculate the change of energy
	  inold = StrengthToGroup(nlist[target],glist[oldg]);
	  innew = StrengthToGroup(nlist[target],glist[newg]);
	  nlink = NodeStrength(nlist[target]);

	  dE = 0.0;

	  dE -= (double)(2 * glist[oldg]->inlinksW) /
		(double)totallinks -
		(double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) *
		(double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) /
		((double)totallinks * (double)totallinks);

	  dE -= (double)(2 * glist[newg]->inlinksW) /
		(double)totallinks -
		(double)(glist[newg]->totlinksW+glist[newg]->inlinksW) *
		(double)(glist[newg]->totlinksW+glist[newg]->inlinksW) /
		((double)totallinks * (double)totallinks);

	  dE += (double)(2*glist[oldg]->inlinksW - 2*inold) /
		(double)totallinks -
		(double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW -
			 nlink ) *
		(double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW -
			 nlink ) /
		((double)totallinks * (double)totallinks);

	  dE += (double)(2*glist[newg]->inlinksW + 2*innew) /
		(double)totallinks -
		(double)(glist[newg]->totlinksW + glist[newg]->inlinksW +
			 nlink ) *
		(double)(glist[newg]->totlinksW + glist[newg]->inlinksW +
			 nlink ) /
		((double)totallinks * (double)totallinks);

	  // Accept the change according to the Boltzman factor
	  if( (dE >= 0.0) || (gsl_rng_uniform(gen) < exp(dE/T)) ){
		MoveNode(nlist[target],glist[oldg],glist[newg]);
		energy += dE;
	  }
	}

	T = T * Ts;
	  } // End of temperature loop
	} // End if totallinks > 0
  }

  RemoveGraph(net);
  return split;
}

// merge = 0 => No group merging
struct group *SACommunityIdentWeight(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int merge, gsl_rng *gen)
{
  int i;
  struct group *part = NULL;
  struct group *split = NULL, *g = NULL;
  struct group **glist;
  struct node_gra **nlist;
  struct node_gra *p;
  struct node_lis *nod;
  int target,empty;
  int newg,oldg;
  int nnod;
  double totallinks = 0.0;
  double innew,inold,nlink;
  double energy = 0.0, dE;
  double T;
  int g1,g2;
  double energyant = 0.0;
  int count = 0, limit = 25; // to stop the search if the energy
						  // does not change
  int cicle1,cicle2;
  int *trans;


  // Create the groups and assign each node to one group
  nnod = CountNodes(net);
  part = CreateHeaderGroup();
  p = net->next;
  ResetNetGroup(net); // All nodes reset to group -1

  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  glist = (struct group **) calloc(nnod, sizeof(struct group *));
  trans = (int *) calloc(nnod,sizeof(int)); 
	
  nlist[0] = p;
  trans[p->num] = 0;
  glist[0] = CreateGroup(part,0);
  AddNodeToGroup(glist[0],p);
  totallinks += NodeStrength(p);

  for( i=1; i<nnod; i++ ) {
	p = p->next;

	nlist[i] = p;
	trans[p->num] = i;
	glist[i] = CreateGroup(glist[i-1],i);
	AddNodeToGroup(glist[i],p);
	totallinks += NodeStrength(p);
  }

  // Number of iterations at each temperature
  if (fac*(double)(nnod*nnod) < 10)
	cicle1 = 10;
  else
	cicle1 = floor(fac*(double)(nnod*nnod));

  if (fac*(double)nnod < 2)
	cicle2 = 2;
  else
	cicle2 = floor(fac*(double)nnod);

  // Do the simulated annealing
  T = Ti;
  energy = ModularityWeight(part);

  while( T > Tf && count < limit){

	if (fabs(energy - energyant) / fabs(energy) < 1.0e-12)
	  count++;
	else{
	  energyant = energy;
	  count = 0;
	}


	if (merge == 1) {

	  for ( i=0; i < cicle2; i++ ){

	//////////////////////////////////////////////////////
	// Propose a pair of merge/split collective changes //
	//////////////////////////////////////////////////////

	// Merge /////////////////////////////////////////////
	target = floor(gsl_rng_uniform(gen) * nnod);
	g1 = nlist[target]->inGroup;

	if(glist[g1]->size < nnod){

	  do{
		target = floor(gsl_rng_uniform(gen) * nnod);
		g2 = nlist[target]->inGroup;
	  }while( g1 == g2 );

	  // Calculate the change of energy
	  nlink = NG2GLinksWeight(glist[g1],glist[g2]);

	  dE = 0.0;

	  dE -= (double)(2*glist[g1]->inlinksW) / (double)totallinks -
		(double)((glist[g1]->totlinksW + glist[g1]->inlinksW) *
			 (glist[g1]->totlinksW + glist[g1]->inlinksW) ) /
		(double)(totallinks * totallinks );

	  dE -= (double)(2*glist[g2]->inlinksW) / (double)totallinks -
		(double)((glist[g2]->totlinksW + glist[g2]->inlinksW) *
			 (glist[g2]->totlinksW + glist[g2]->inlinksW) ) /
		(double)(totallinks * totallinks );

	  dE += 2.0*(double)(glist[g1]->inlinksW +
				 glist[g2]->inlinksW+nlink)
		/ (double)totallinks -
		(double)(glist[g1]->totlinksW + glist[g1]->inlinksW +
			 glist[g2]->totlinksW + glist[g2]->inlinksW ) *
		(double)(glist[g1]->totlinksW + glist[g1]->inlinksW +
			 glist[g2]->totlinksW + glist[g2]->inlinksW ) /
		(double)(totallinks*totallinks);

	  // Accept the change according to Metroppolis
	  if( (dE >= 0.0) || ( gsl_rng_uniform(gen) < exp(dE/T) ) ){
		MergeGroups(glist[g1],glist[g2]);
		energy += dE;
	  }
	}

	// Split /////////////////////////////////////////////
	target = floor(gsl_rng_uniform(gen) * nnod); // target node
	target = nlist[target]->inGroup;    // target group

	// Look for an empty group
	g = part;
	empty = -1;
	while((g->next != NULL) && (empty < 0)){
	  g = g->next;
	  if (g->size == 0){
		empty = g->label;
	  }
	}

	if (empty >= 0 ){ // if there are no empty groups, do nothing
/*	  split = BestNetworkSplitWeight(glist[target],gen); */
/*	  split = ThermalNetworkSplitWeight(glist[target],Ti,T,gen); */
	  split = ThermalPercNetworkSplitWeight(glist[target],
						Ti, T, gen);

	  // Split the group
	  nod = (split->next)->nodeList;
	  while ( nod->next != NULL ){
		nod = nod->next;
		MoveNode(nlist[trans[nod->node]],
			 glist[target],
			 glist[empty]);
	  }
	  RemovePartition(split);
	  split = NULL;

	  // Try to re-merge the two groups
	  // Calculate the change of energy
	  nlink = NG2GLinksWeight(glist[target],glist[empty]);

	  dE = 0.0;

	  dE -= (double)(2*glist[target]->inlinksW) /
		(double)totallinks -
		(double)((glist[target]->totlinksW +
			  glist[target]->inlinksW) *
			 (glist[target]->totlinksW +
			  glist[target]->inlinksW) ) /
		(double)(totallinks * totallinks );

	  dE -= (double)(2*glist[empty]->inlinksW) /
		(double)totallinks -
		(double)((glist[empty]->totlinksW +
			  glist[empty]->inlinksW) *
			 (glist[empty]->totlinksW +
			  glist[empty]->inlinksW) ) /
		(double)(totallinks * totallinks );

	  dE += 2.0*(double)(glist[target]->inlinksW +
				 glist[empty]->inlinksW+nlink)
		/ (double)totallinks -
		(double)(glist[target]->totlinksW +
			 glist[target]->inlinksW +
			 glist[empty]->totlinksW +
			 glist[empty]->inlinksW) *
		(double)(glist[target]->totlinksW +
			 glist[target]->inlinksW +
			 glist[empty]->totlinksW +
			 glist[empty]->inlinksW) /
		(double)(totallinks*totallinks);

	  // Accept the change according to "inverse" Metroppolis.
	  // Inverse means that the algor is applied to the split
	  // and NOT to the merge!
	  if( (dE >= 0.0) && ( gsl_rng_uniform(gen) > exp(-dE/T) ) ){
		MergeGroups(glist[target],glist[empty]);
	  }
	  else{
		energy -= dE;
	  }

	} // End of if empty

	  } // End of cicle2 loop
	} // End of if merge == 1

	for ( i=0; i < cicle1; i++ ){

	  ///////////////////////////////
	  // Propose an individual change
	  ///////////////////////////////
	  target = floor(gsl_rng_uniform(gen) * (double)nnod);
	  oldg = nlist[target]->inGroup;
	  do{
	newg = floor(gsl_rng_uniform(gen) * (double)nnod);
	  }while(newg == oldg);

	  // Calculate the change of energy
	  inold = StrengthToGroup(nlist[target],glist[oldg]);
	  innew = StrengthToGroup(nlist[target],glist[newg]);
	  nlink = NodeStrength(nlist[target]);

	  dE = 0.0;

	  dE -= (double)(2 * glist[oldg]->inlinksW) /
	(double)totallinks -
	(double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) *
	(double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) /
	((double)totallinks * (double)totallinks);

	  dE -= (double)(2 * glist[newg]->inlinksW) /
	(double)totallinks -
	(double)(glist[newg]->totlinksW+glist[newg]->inlinksW) *
	(double)(glist[newg]->totlinksW+glist[newg]->inlinksW) /
	((double)totallinks * (double)totallinks);

	  dE += (double)(2*glist[oldg]->inlinksW - 2*inold) /
	(double)totallinks -
	(double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW -
		 nlink ) *
	(double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW -
		 nlink ) /
	((double)totallinks * (double)totallinks);

	  dE += (double)(2*glist[newg]->inlinksW + 2*innew) /
	(double)totallinks -
	(double)(glist[newg]->totlinksW + glist[newg]->inlinksW +
		 nlink ) *
	(double)(glist[newg]->totlinksW + glist[newg]->inlinksW +
		 nlink ) /
	((double)totallinks * (double)totallinks);

	  // Accept the change according to Metroppolis
	  if( (dE >= 0.0) || ( gsl_rng_uniform(gen) < exp(dE/T) ) ){
	MoveNode(nlist[target],glist[oldg],glist[newg]);
	energy += dE;
	  }
	}

	T = T * Ts;
  }

/*   printf("energy = %lf\n",energy); */

  return CompressPart(part);
}

/*
  ---------------------------------------------------------------------
  Splits one group into two new groups (thermalized with the outer
  SA), checking first if the network contains disconnected clusters.
  ---------------------------------------------------------------------
*/
void
SAGroupSplitBipart(struct group *target_g, struct group *empty_g,
		   double Ti, double Tf, double Ts,
		   double cluster_prob,
		   double **cmat, double msfac,
		   gsl_rng *gen)
{
  struct group *glist[2], *g = NULL, *split = NULL;
  struct node_gra **nlist;
  struct node_lis *p = NULL;
  int nnod = 0;
  int i;
  int n1, n2, t1, t2;
  int target, oldg, newg;
  double dE = 0.0;
  double T;
  double dice;
  struct node_gra *net;
  struct binet *binet;
  int S1 = target_g->size;
  int cluster_sw;
  void *dict=NULL;
  double energy, energyant;
  int count = 0, limit = 5;

  /* Initialize */
  glist[0] = target_g;
  glist[1] = empty_g;

  /* Are we going to use clusters? */
  if (gsl_rng_uniform(gen) < cluster_prob)
	cluster_sw = 1;

  /*
	If cluster_sw, check if the network is disconnected
  */
  if (cluster_sw == 1) {
	/* Build a network from the nodes in the target group and find
	   disconected clusters  */
	binet = CreateBipart();
	binet->net1 = BuildNetFromGroup(target_g);
	binet->net2 = CreateHeaderGraph();
	net = ProjectBipart(binet); /* This trick to generate the projection
				  works, although it is not elegant */
	split = ClustersPartition(net);
  }

  /*
	If cluster_sw==1 and the network is disconnected, then use the
	clusters
  */
  if (cluster_sw == 1 && NGroups(split) > 1) {
	/* Build a dictionary for fast access to nodes */
	dict = MakeLabelDict(binet->net1);

	/* Move the nodes in some (half) of the clusters to empty module */
	g = split;
	while ((g = g->next) != NULL) {
	  if (gsl_rng_uniform(gen) < 0.500000) { // With prob=0.5 move to empty
	p = g->nodeList;
	while ((p = p->next) != NULL) {
	  MoveNode(GetNodeDict(p->nodeLabel, dict), target_g, empty_g);
	}
	  }
	}
  }

  /*
	Else if the network is connected or cluster_sw... go SA!
  */
  else {
	/* Allocate memory for a list of nodes for faster access */
	nlist = (struct node_gra **) calloc(S1, sizeof(struct node_gra *));
	nnod = 0;

	/* Randomly assign the nodes to the groups */
	p = target_g->nodeList;
	while (p->next != NULL) {
	  nlist[nnod++] = p->next->ref;
	  dice = gsl_rng_uniform(gen);
	  if (dice < 0.500000) {
	MoveNode(p->next->ref, target_g, empty_g);
	  }
	  else {
	p = p->next;
	  }
	}

	/* Do SA to "optimize" the splitting */
	T = Ti;
	energy = energyant = 0.0;
	while ((T >= Tf) && (count < limit)) {

	  /* Do nnod moves */
	  for (i=0; i<nnod; i++) {

	/* Determine target node */
	target = floor(gsl_rng_uniform(gen) * (double)nnod);
	if (nlist[target]->inGroup == target_g->label)
	  oldg = 0;
	else
	  oldg = 1;
	newg = 1 - oldg;

	/* Calculate the change of energy */
	dE = 0.0;
	n1 = nlist[target]->num;
	t1 = CountLinks(nlist[target]);
	/* 1-Old group */
	p = glist[oldg]->nodeList;
	while ((p = p->next) != NULL) {
	  n2 = p->node;
	  if (n2 != n1) {
		t2 = CountLinks(p->ref);
		dE -= 2. * (cmat[n1][n2] - t1 * t2 * msfac);
	  }
	}
	/* 2-New group */
	p = glist[newg]->nodeList;
	while ((p = p->next) != NULL) {
	  n2 = p->node;
	  t2 = CountLinks(p->ref);
	  dE += 2. * (cmat[n1][n2] - t1 * t2 * msfac);
	}

	/* Accept the change according to the Boltzman factor */
	if (gsl_rng_uniform(gen) < exp(dE/T)) {
	  MoveNode(nlist[target], glist[oldg], glist[newg]);
	  energy += dE;
	}
	  }

	  /* Update the no-change counter */
	  if (fabs(energy - energyant) / fabs(energyant) < EPSILON_MOD_B ||
	  fabs(energyant) < EPSILON_MOD_B)
	count++;
	  else
	count = 0;

	  /* Update the last energy */
	  energyant = energy;

	  /* Update the temperature */
	  T = T * Ts;

	} /* End of temperature loop */

  } /* End of else (network is connected) */

  if (cluster_sw == 1) {
	RemovePartition(split);
	RemoveGraph(net);
	RemoveBipart(binet);
	FreeLabelDict(dict);
  }
  else {
	free(nlist);
  }

  /* Done */
  return;
}

/*
  ---------------------------------------------------------------------
  Splits one group into two new groups (thermalized with the outer
  SA), checking first if the network contains disconnected clusters.
  For Weighted Bipartite Networks.
  ---------------------------------------------------------------------
*/
void
SAGroupSplitBipartWeighted(struct group *target_g, struct group *empty_g,
			   double Ti, double Tf, double Ts,
			   double cluster_prob,
			   double *strength,
			   double **swwmat, double Wafac,
			   gsl_rng *gen)
{
  struct group *glist[2], *g = NULL, *split = NULL;
  struct node_gra **nlist;
  struct node_lis *p = NULL;
  int nnod = 0;
  int i;
  int n1, n2;
  int target, oldg, newg;
  double dE = 0.0;
  double T;
  double dice;
  struct node_gra *net;
  struct binet *binet;
  int SZ1 = target_g->size;
  int cluster_sw;
  void *dict=NULL;
  double energy, energyant;
  int count = 0, limit = 5;

  /* Initialize */
  glist[0] = target_g;
  glist[1] = empty_g;

  /* Are we going to use clusters? */
  if (gsl_rng_uniform(gen) < cluster_prob)
	cluster_sw = 1;

  /*
	If cluster_sw, check if the network is disconnected
  */
  if (cluster_sw == 1) {
	/* Build a network from the nodes in the target group and find disconected clusters  */
	binet = CreateBipart();
	binet->net1 = BuildNetFromGroup(target_g);
	binet->net2 = CreateHeaderGraph();
	net = ProjectBipart(binet);
	/* This trick to generate the projection works, although it is not elegant */
	/*
	  Note: This projection uses NCommonLinksBipart to determine
	  weights in the projected network, and the original, bipartite
	  weights are disregarded. But that's ok, since the projected
	  network is just used to find clusters.
	*/
	split = ClustersPartition(net);
  }

  /*
	If cluster_sw==1 and the network is disconnected, then use the clusters
  */
  if (cluster_sw == 1 && NGroups(split) > 1) {
	/* Build a dictionary for fast access to nodes */
	dict = MakeLabelDict(binet->net1);

	/* Move the nodes in some (that is, approximately half) of the clusters to the empty group */
	g = split;
	while ((g = g->next) != NULL) {
	  if (gsl_rng_uniform(gen) < 0.500000) { // With prob=0.5 move to empty
		p = g->nodeList;
		while ((p = p->next) != NULL) {
		  MoveNodeFast(GetNodeDict(p->nodeLabel, dict), target_g, empty_g);
		}
	  }
	}
  }

  /*
	Else if the network is connected or cluster_sw... go SA!
  */
  else {
	/* Allocate memory for a list of nodes for faster access */
	nlist = (struct node_gra **) calloc(SZ1, sizeof(struct node_gra *));
	nnod = 0;

	/* Randomly assign the nodes to the groups */
	p = target_g->nodeList;
	while (p->next != NULL) {
	  nlist[nnod++] = p->next->ref;
	  dice = gsl_rng_uniform(gen);
	  if (dice < 0.500000) {
		MoveNodeFast(p->next->ref, target_g, empty_g);
	  }
	  else {
		p = p->next;
	  }
	}

	/* Do SA to "optimize" the splitting */
	T = Ti;
	energy = energyant = 0.0;
	while ((T >= Tf) && (count < limit)) {

	  /* Do nnod moves */
	  for (i=0; i<nnod; i++) {
		/* Determine target node */
		target = floor(gsl_rng_uniform(gen) * (double)nnod);
		if (nlist[target]->inGroup == target_g->label)
		  oldg = 0;
		else
		  oldg = 1;
		newg = 1 - oldg;

		  /* Calculate the change of energy */
		dE = 0.0;
		n1 = nlist[target]->num;

		/* 1-Old group */
		p = glist[oldg]->nodeList;
		while ((p = p->next) != NULL) {
			n2 = p->node;
		  if (n2 != n1) {
			dE -= swwmat[n1][n2];
		  }
		}

		/* 2-New group */
		p = glist[newg]->nodeList;
		while ((p = p->next) != NULL) {
		  n2 = p->node;
		  dE += swwmat[n1][n2];
		}

		  /* Accept the change according to the Boltzman factor */
		if (gsl_rng_uniform(gen) < exp(dE/T)) {
		  MoveNodeFast(nlist[target], glist[oldg], glist[newg]);
		  energy += dE;
		}
	  }

	  /* Update the no-change counter */
	  if (fabs(energy - energyant) / fabs(energyant) < EPSILON_MOD_B || fabs(energyant) < EPSILON_MOD_B)
		  count++;
	  else
		count = 0;

	  /* Update the last energy */
	  energyant = energy;

	  /* Update the temperature */
	  T = T * Ts;

	} /* End of temperature while loop */

  } /* End of else (network is connected) */

  if (cluster_sw == 1) {
	RemovePartition(split);
	RemoveGraph(net);
	RemoveBipart(binet);
	FreeLabelDict(dict);
  }
  else {
	free(nlist);
  }

  /* Done */
  return;
}

/*
  ---------------------------------------------------------------------
  Identify modules in the net1 network of a bipartite network using
  simulated annealing.

  Ti: Initial temperature for the SA

  Tf: Final temperature for the SA

  Ts: Cooling factor

  fac: Iteration factor

  merge: implement collective moves (merge=1) or not (merge=0)

  prob: only used if merge==1. Use percolation split with probability
  prob (prob>0) or do not use percolation (prob<0)
  ---------------------------------------------------------------------
*/
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
  int i;
  struct node_gra *net1 = binet->net1;
  struct node_gra *net2 = binet->net2;
  double sms = 0.0, sms2 = 0.0, msfac;
  struct node_gra *p, *p2;
  int nnod;
  struct group *part = NULL, *g = NULL;
  struct node_gra **nlist;
  struct group **glist = NULL, *lastg;
  int cicle1, cicle2;
  int count = 0, limit = 25;
  double energy, energyant=0.0, dE;
  double T;
  int target, oldg, newg;
  double **cmat;
  int g1, g2, empty;
  struct node_lis *nod, *nod2;
  double t1, t2;
  struct group *best_part = NULL;
  double best_E = -100.0;
  double cluster_prob = 0.500000;
  int dice;
  void *nodeDict;
  int *nlink=NULL;
  FILE *outf;

  /*
	Preliminaries: Initialize, allocate memory, and place nodes in
	initial groups
	-------------------------------------------------------------------
  */
  /* Create the groups and assign each node to one group */
  nnod = CountNodes(net1);
  ResetNetGroup(net1);
  part = CreateHeaderGroup();
  p = net1;

  /* Create a node dictionary for fast access to nodes by label */
  nodeDict = MakeLabelDict(net1);

  /* Allocate memory for the node list */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));

  /* Create the groups and assign nodes to the initial group according
	 to initial_sw. Additionally, map nodes and groups to lists for
	 faster access. */
  switch (initial_sw) {

  case 'o':         /*  One node in each group */
	ngroup = nnod;
	glist = (struct group **) calloc(ngroup, sizeof(struct group *));
	while ((p = p->next) != NULL) {
	  glist[p->num] = CreateGroup(part, p->num);
	  nlist[p->num] = p;
	  AddNodeToGroup(glist[p->num], p);
	}
	break;

  case 'r':        /*  Random placement of nodes in groups */
	if (ngroup < 1)
	  ngroup = nnod;
	glist = (struct group **) calloc(ngroup, sizeof(struct group *));
	lastg = part;
	for (i=0; i<ngroup; i++)
	  glist[i] = lastg = CreateGroup(lastg, i);
	while ((p = p->next) != NULL) {
	  nlist[p->num] = p;
	  dice = floor(gsl_rng_uniform(gen)* (double)ngroup);
	  AddNodeToGroup(glist[dice], p);
	}
	break;
  }

  /* Calculate s2ms=(sum m_s)^2 and sms2=sum(m_s^2) */
  p = net2;
  while ((p = p->next) != NULL) {
	sms += (double)CountLinks(p);
	sms2 += (double)(CountLinks(p) * CountLinks(p));
  }
  msfac = 1. / (sms * sms);

  /* Calculate the c matrix of concurrences, and the number of links
	 of each node */
  cmat = allocate_d_mat(nnod, nnod);
  nlink = allocate_i_vec(nnod);
  p = net1;
  while ((p = p->next) != NULL) {
	nlink[p->num] = CountLinks(p);
	p2 = net1;
	while ((p2 = p2->next) != NULL) {
	  cmat[p->num][p2->num] = (double)NCommonLinksBipart(p, p2) /
	(sms2 - sms);
	}
  }

  /*
	Determine the number of iterations at each temperature
	-------------------------------------------------------------------
  */
  if (fac * (double)(nnod * nnod) < 10)
	cicle1 = 10;
  else
	cicle1 = floor(fac * (double)(nnod * nnod));

  if (fac * (double)nnod < 2)
	cicle2 = 2;
  else
	cicle2 = floor(fac * (double)nnod);

  /*
	Do the simulated annealing
	-------------------------------------------------------------------
  */
  /* Determine initial values */
  T = Ti;
  energy = ModularityBipart(binet, part);

  /* Temperature loop */
  while ((T >= Tf) && (count < limit)) {

	/* Output */
	switch (output_sw) {
	case 'n':
	  break;
	case 'b':
	  break;
	case 's':
	  fprintf(stderr, "%g %lf %g\n",1.0/T, energy, T);
	  break;
	case 'm':
	  fprintf(stderr, "%g %lf %g\n",1.0/T, energy, T);
	  break;
	case 'v':
	  fprintf(stderr, "%g %lf %lf %g\n",
		  1.0/T, energy, ModularityBipart(binet, part), T);
	  break;
	case 'd':
	  FPrintPartition(stderr, part, 0);
	  fprintf(stderr, "%g %lf %lf %g\n",
		  1.0/T, energy, ModularityBipart(binet, part), T);
	}


	/*
	  Do cicle1 individual change iterations
	*/
	for (i=0; i<cicle1; i++) {

	  /* Propose an individual change */
	  target = floor(gsl_rng_uniform(gen) * (double)nnod);
	  oldg = nlist[target]->inGroup;
	  do {
	newg = floor(gsl_rng_uniform(gen) * ngroup);
	  } while (newg == oldg);

	  /* Calculate the change of energy */
	  dE = 0.0;
	  t1 = nlink[target];

	  /* Old group contribution */
	  nod = glist[oldg]->nodeList;
	  while ((nod = nod->next) != NULL) {
	t2 = nlink[nod->ref->num];
	dE -= 2. * (cmat[nlist[target]->num][nod->node] -
			t1 * t2 * msfac);
	  }

	  /* New group contribution */
	  nod = glist[newg]->nodeList;
	  while ((nod = nod->next) != NULL) {
	t2 = nlink[nod->ref->num];
	dE += 2. * (cmat[nlist[target]->num][nod->node] -
			t1 * t2 * msfac);
	  }
	  dE += 2. * (cmat[nlist[target]->num][nlist[target]->num] -
		  t1 * t1 * msfac);

	  /* Accept or reject movement according to Metropolis */
	  if (gsl_rng_uniform(gen) < exp(dE/T)) {
	energy += dE;
	MoveNode(nlist[target],glist[oldg],glist[newg]);
	  }
	}

	/*
	  Do cicle2 collective change iterations
	*/
	if (collective_sw == 1) {
	  for (i=0; i<cicle2; i++){

	/* MERGE */
	target = floor(gsl_rng_uniform(gen) * nnod);
	g1 = nlist[target]->inGroup;

	if (glist[g1]->size < nnod) {
	  do {
		target = floor(gsl_rng_uniform(gen) * nnod);
		g2 = nlist[target]->inGroup;
	  } while (g1 == g2);

	  /* Calculate dE */
	  dE = 0.0;
	  nod = glist[g1]->nodeList;
	  while ((nod = nod->next) != NULL) {
		nod2 = glist[g2]->nodeList;
		while ((nod2 = nod2->next) != NULL) {
		  t1 = nlink[nod->ref->num];
		  t2 = nlink[nod2->ref->num];
		  dE += 2. * (cmat[nod->node][nod2->node] -
			  t1 * t2 * msfac);
		}
	  }

	  /* Accept or reject change */
	  if (gsl_rng_uniform(gen) < exp(dE/T)) {
		MergeGroups(glist[g1], glist[g2]);
		energy += dE;
	  }
	} /* End of merge move */

	/* SPLIT */
	/* Look for an empty group */
	g = part;
	empty = -1;
	while (((g = g->next) != NULL) && (empty < 0)) {
	  if (g->size == 0) {
		empty = g->label;
		break;
	  }
	}

	if (empty >= 0 ) { /* if there are no empty groups, do nothing */
	  /* Select group to split */
	  do {
		target = floor(gsl_rng_uniform(gen) * (double)nnod); /* node */
		target = nlist[target]->inGroup;    /* target group */
	  } while (glist[target]->size == 1);

	  /* Split the group */
	  SAGroupSplitBipart(glist[target], glist[empty],
				 Ti, T, 0.95,
				 cluster_prob,
				 cmat, msfac, gen);

	  /* Calculate dE for remerging the groups */
	  dE = 0.0;
	  nod = glist[target]->nodeList;
	  while ((nod = nod->next) != NULL) {
		nod2 = glist[empty]->nodeList;
		while ((nod2 = nod2->next) != NULL) {
		  t1 = nlink[nod->ref->num];
		  t2 = nlink[nod2->ref->num];
		  dE += 2. * (cmat[nod->node][nod2->node] -
			  t1 * t2 * msfac);
		}
	  }

	  /* Accept the change according to "inverse" Metroppolis.
		 Inverse means that the algor is applied to the split and
		 NOT to the merge! */
	  if ((dE > EPSILON_MOD_B) && (gsl_rng_uniform(gen) > exp(-dE/T))) {
		/* Undo the split */
		MergeGroups(glist[target],glist[empty]);
	  }
	  else{
		/* Update energy */
		energy -= dE;
	  }
	} /* End of split move */
	  } /* End of cicle2 loop */
	} /* End of 'if collective_sw==1' loop */

	/* Update the no-change counter */
	if (((T < Ti / 1000.) || (Ti < EPSILON_MOD_B)) &&
	(fabs(energy - energyant) / fabs(energyant) < EPSILON_MOD_B ||
	fabs(energyant) < EPSILON_MOD_B)) {
	  count++;

	  /* If the SA is ready to stop (count==limit) but the current
	 partition is not the best one so far, replace the current
	 partition by the best one and continue from there. */
	  if ((count == limit) && (energy + EPSILON_MOD_B < best_E)) {
	switch (output_sw) {
	case 'n':
	  break;
	case 'b':
	  break;
	default:
	  fprintf(stderr, "# Resetting partition\n");
	  break;
	}

	/* Remap the partition to the best partition */
	RemovePartition(part);
	g = part = CopyPartition(best_part);
	while ((g = g->next) != NULL)
	  glist[g->label] = g;
	MapPartToNet(part, binet->net1);

	/* Reset energy and counter */
	energy = best_E;
	count = 0;
	  }
	}

	else {
	  count = 0;
	}
	/* Update the last energy */
	energyant = energy;

	/* Compare the current partition to the best partition so far and
	   save the current if it is better than the best so far. */
	if ( energy > best_E ) {
	  if ( best_part != NULL )
	RemovePartition(best_part);
	  best_part = CopyPartition(part);
	  MapPartToNet(part, binet->net1); /* MUST DO this after copying a
					  part! */
	  best_E = energy;
	}

	/* Save the partition to a file if necessary */
	switch (output_sw) {
	case 'b':
	  outf = fopen("part.tmp", "w");
	  FPrintPartition(outf, best_part, 1);
	  fclose(outf);
	case 's':
	  outf = fopen("part.tmp", "w");
	  FPrintPartition(outf, best_part, 1);
	  fclose(outf);
	default:
	  break;
	}

	/* Update the temperature */
	T = T * Ts;

  } /* End of simulated annealing */

  /* Free memory */
  free_d_mat(cmat, nnod);
  free_i_vec(nlink);
  RemovePartition(best_part);
  FreeLabelDict(nodeDict);
  free(glist);
  free(nlist);

  /* Done */
  return CompressPart(part);
}

/*
  ---------------------------------------------------------------------
  Identify modules in the net1 network of a WEIGHTED bipartite network
  using simulated annealing.

  Ti: Initial temperature for the SA

  Tf: Final temperature for the SA

  Ts: Cooling factor

  fac: Iteration factor

  merge: implement collective moves (merge=1) or not (merge=0)

  prob: only used if merge==1. Use percolation split with probability
  prob (prob>0) or do not use percolation (prob<0)
  ---------------------------------------------------------------------
*/
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
  int i;
  struct node_gra *net1 = binet->net1;
  struct node_gra *net2 = binet->net2;
  double sWa = 0.0, s2Wa = 0.0, sWa2 = 0.0, Wafac, s;
  struct node_gra *p, *p2;
  int nnod;
  struct group *part = NULL, *g = NULL;
  struct node_gra **nlist;
  struct group **glist = NULL, *lastg;
  int cicle1, cicle2;
  int count = 0, limit = 25;
  double energy, energyant = 0.0, dE;
  double T;
  int target, oldg, newg;
  char* accepted;
  double **swwmat;
  int g1, g2, empty;
  struct node_lis *nod, *nod2;
  double s1, s2, sww12;
  struct group *best_part = NULL;
  double best_E = -100.0;
  double cluster_prob = 0.500000;
  int dice;
  void *nodeDict;
  double *strength=NULL;
  FILE *outf;

  /*
	Preliminaries: Initialize, allocate memory, and place nodes in
	initial groups
	-------------------------------------------------------------------
  */
  /* Create the groups and assign each node to one group */
  nnod = CountNodes(net1);
  ResetNetGroup(net1);
  part = CreateHeaderGroup();
  p = net1;

  /* Create a node dictionary for fast access to nodes by label */
  nodeDict = MakeLabelDict(net1);

  /* Allocate memory for the node list */
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));

  /* Create the groups and assign nodes to the initial group according
	 to initial_sw. Additionally, map nodes and groups to lists for
	 faster access. */
  switch (initial_sw) {

  case 'o':         /*  One node in each group */
	ngroup = nnod;
	glist = (struct group **) calloc(ngroup, sizeof(struct group *));
	while ((p = p->next) != NULL) {
	  glist[p->num] = CreateGroup(part, p->num);
	  nlist[p->num] = p;
	  AddNodeToGroupFast(glist[p->num], p);
	}
	break;

  case 'r':        /*  Random placement of nodes in groups */
	if (ngroup < 1)
	  ngroup = nnod;
	glist = (struct group **) calloc(ngroup, sizeof(struct group *));
	lastg = part;
	for (i=0; i<ngroup; i++)
	  glist[i] = lastg = CreateGroup(lastg, i);
	while ((p = p->next) != NULL) {
	  nlist[p->num] = p;
	  dice = floor(gsl_rng_uniform(gen)* (double)ngroup);
	  AddNodeToGroupFast(glist[dice], p);
	}
	break;
  }


  /* Calculate s2Wa=(sum W_a)^2 and sWa2=sum(W_a^2) */
  p = net2;
  while ((p = p->next) != NULL) {
	s = (double)NodeStrengthFast(p);
	sWa += s;
	sWa2 += s * s;
  }
  s2Wa = sWa * sWa;
  Wafac = 1. / s2Wa;

  /*
	 Calculate the matrix of all pairwise contributions to modularity
  */
  swwmat = allocate_d_mat(nnod, nnod);
  p = net1;
  while ((p = p->next) != NULL) {
	s1 = (double)NodeStrengthFast(p);
	p2 = net1;
	while ((p2 = p2->next) != NULL) {
	  s2 = (double)NodeStrengthFast(p2);
	  if(swwmat[p2->num][p->num] != 0)
		swwmat[p->num][p2->num] = swwmat[p2->num][p->num];
	  else{
		sww12 = (double)SumProductsOfCommonWeightsBipart(p, p2);
		swwmat[p->num][p2->num] = 2. * (sww12 / (sWa2 - sWa) - s1 * s2 / s2Wa);
	  }
	}
  }


  /*
	Determine the number of iterations at each temperature
	-------------------------------------------------------------------
  */
  if (fac * (double)(nnod * nnod) < 10)
	cicle1 = 10;
  else
	cicle1 = floor(fac * (double)(nnod * nnod));

  if (fac * (double)nnod < 2)
	cicle2 = 2;
  else
	cicle2 = floor(fac * (double)nnod);


  /*
	Do the simulated annealing
	-------------------------------------------------------------------
  */
  /* Determine initial values */
  T = Ti;
  energy = ModularityBipartWeightedFast(binet, part, swwmat);

  /* Temperature loop */
  while ((T >= Tf) && (count < limit)) {

	/* Output */
	switch (output_sw) {
	case 'n':
	  break;
	case 'b':
	  break;
	case 's':
	  fprintf(stderr, "%g %lf %g\n",1.0/T, energy, T);
	  break;
	case 'm':
	  fprintf(stderr, "%g %lf %g\n",1.0/T, energy, T);
	  break;
	case 'v':
	  fprintf(stderr, "%g %lf %lf %g\n", 1.0/T, energy, ModularityBipartWeightedFast(binet, part, swwmat), T);
	  break;
	case 'd':
	  FPrintPartition(stderr, part, 0);
	  fprintf(stderr, "%g %lf %lf %g\n", 1.0/T, energy, ModularityBipartWeightedFast(binet, part, swwmat), T);
	}

	/*
	  Do cicle1 individual change iterations
	*/
	for (i=0; i<cicle1; i++) {

	  accepted = "REJECTED";

	  /* Propose an individual change */
	  target = floor(gsl_rng_uniform(gen) * (double)nnod);
	  oldg = nlist[target]->inGroup;
	  do {
		newg = floor(gsl_rng_uniform(gen) * ngroup);
	  } while (newg == oldg);

	  /* Calculate the change of energy */
	  dE = 0.0;

	  /* Old group contribution */
	  nod = glist[oldg]->nodeList;
	  while ((nod = nod->next) != NULL) {
		if(nlist[target]->num != nod->node){
		  dE -= swwmat[nlist[target]->num][nod->node];
		}
	  }

	  /* New group contribution */
	  nod = glist[newg]->nodeList;
	  while ((nod = nod->next) != NULL) {
		dE += swwmat[nlist[target]->num][nod->node];
	  }

	  /* Accept or reject movement according to Metropolis */
	  if (gsl_rng_uniform(gen) < exp(dE/T)) {
		accepted = "ACCEPTED";
		energy += dE;
		MoveNodeFast(nlist[target],glist[oldg],glist[newg]);
	  }

	  switch(output_sw) {
	  case 'd':
		if (dE < 0 && strcmp(accepted,"ACCEPTED")==0 && (T < Ti/1.0e6)) {
		  fprintf(stderr, "Cicle 1 move %s: %i move %i -> %i dE=%g prob=%g T=%g\n", accepted, target, oldg, newg, dE, exp(dE/T),T);
		}
		  accepted = "REJECTED";
	  }
	}

	/*
	  Do cicle2 collective change iterations
	*/
	if (collective_sw == 1) {
	  for (i=0; i<cicle2; i++){

		/* MERGE */
		accepted = "REJECTED";

		target = floor(gsl_rng_uniform(gen) * nnod);
		g1 = nlist[target]->inGroup;

		/* Can only merge groups when there isn't a single megagroup */
		if (glist[g1]->size < nnod) {
		  do {
			target = floor(gsl_rng_uniform(gen) * nnod);
			g2 = nlist[target]->inGroup;
		  } while (g1 == g2);

			/* Calculate dE */
		  dE = 0.0;
		  nod = glist[g1]->nodeList;
		  while ((nod = nod->next) != NULL) {
			nod2 = glist[g2]->nodeList;
			while ((nod2 = nod2->next) != NULL) {
			  dE += swwmat[nod->node][nod2->node];
			}
		  }

			/* Accept or reject change */
		  if (gsl_rng_uniform(gen) < exp(dE/T)) {
			accepted = "ACCEPTED";
			MergeGroupsFast(glist[g1], glist[g2]);
			energy += dE;
		  }

			switch(output_sw) {
		  case 'd':
			if (dE < 0 && strcmp(accepted,"ACCEPTED")==0 && (T < Ti/1.0e6)) {
				fprintf(stderr, "Cicle 2 Merge %s: %i (sz %i) with %i (sz %i) dE=%g prob=%g T=%g\n",
					  accepted, g1, glist[g1]->size, g2, glist[g2]->size, dE, exp(dE/T),T);
			  }
			  accepted = "REJECTED";
			}
		  } /* End of merge move */

		  /* SPLIT */
		  accepted = "REJECTED";

		  /* Look for an empty group to split the target group into*/
		  g = part;
		  empty = -1;
		  while (((g = g->next) != NULL) && (empty < 0)) {
			if (g->size == 0) {
			empty = g->label;
			break;
			}
		  }

		  if (empty >= 0 ) { /* if there are no empty groups, do nothing */
			/* Select group to split */
			do {
			target = floor(gsl_rng_uniform(gen) * (double)nnod); /* node */
			target = nlist[target]->inGroup;    /* target group */
			} while (glist[target]->size == 1);

			/* Split the group */
			SAGroupSplitBipartWeighted(glist[target], glist[empty],
											 Ti, T, 0.95,
									 cluster_prob,
									 strength,
									 swwmat, Wafac, gen);

			/* Calculate dE for remerging the groups */
		  dE = 0.0;
		  nod = glist[target]->nodeList;
		  while ((nod = nod->next) != NULL) {
			nod2 = glist[empty]->nodeList;
			while ((nod2 = nod2->next) != NULL) {
			  dE += swwmat[nod->node][nod2->node];
			  }
			}

			/*
			 Accept the change according to "inverse" Metroppolis.
			   Inverse means that the algor is applied to the split and
			   NOT to the merge!
		  */
			if (gsl_rng_uniform(gen) > exp(-dE/T)) {
			  /* Undo the split */
			  MergeGroupsFast(glist[target],glist[empty]);
			}
			else{
			  accepted = "ACCEPTED";
			  /* Update energy */
			  energy -= dE;
			}

			switch(output_sw) {
			case 'd':
			  if (-dE < 0 && strcmp(accepted,"ACCEPTED")==0  && (T < Ti/1.0e6)) {
				fprintf(stderr, "Cicle 2 Split %s: %i -> %i (sz %i) and %i (sz %i) dE=%g dE>-en*Epsilon/10=%i prob=%g T=%g\n",
						  accepted, target, target, glist[target]->size, empty, glist[empty]->size, -dE, -dE>-fabs(energyant)*EPSILON_MOD_B/10, exp(-dE/T), T);
			  }
			  accepted = "REJECTED";
			}
		  } /* End of split move */
	  } /* End of cicle2 loop */
	} /* End of 'if collective_sw==1' loop */

	/* Update the no-change counter */
	// condition (T < Ti / 1000.) && removed (stpdescent)
	if ((fabs(energy - energyant) / fabs(energyant) < EPSILON_MOD_B || fabs(energyant) < EPSILON_MOD_B)) {
	  count++;


	  /*
		 If the SA is ready to stop (count==limit) but the current
		 partition is not the best one so far, replace the current
		 partition by the best one and continue from there.
	  */
	  if ((count == limit) && (output_sw == 'd')) {
		fprintf(stderr, "# Limit reached.\n");
	  }

	  if ((count == limit) && (energy + EPSILON_MOD_B < best_E)) {
		switch (output_sw) {
		case 'n':
		  break;
		case 'b':
		  break;
		default:
		  fprintf(stderr, "# Resetting partition\n");
			break;
		}

		  /* Remap the partition to the best partition */
		  RemovePartition(part);
		  g = part = CopyPartition(best_part);
		  while ((g = g->next) != NULL)
			glist[g->label] = g;
		  MapPartToNetFast(part, binet->net1);

		  /* Reset energy and counter */
		  energy = best_E;
		  count = 0;
	  }
	}
	else {
	  count = 0;
	}

	if ((T < Ti / 1.0e6) && (output_sw == 'd')) {
	  fprintf(stderr, "En Change %%: %g (Last En: %g), count: %i",
				(fabs(energy - energyant)/fabs(energyant)), fabs(energyant), count);
	}

	/* Update the last energy */
	energyant = energy;

	/* Compare the current partition to the best partition so far and
	   save the current if it is better than the best so far. */
	if ( energy > best_E ) {
	  if ( best_part != NULL )
		  RemovePartition(best_part);
	  best_part = CopyPartition(part);
	  MapPartToNetFast(part, binet->net1); /* MUST DO this after copying a part! */
	  best_E = energy;
	}

	/* Save the partition to a file if necessary */
	switch (output_sw) {
	case 'b':
	  outf = fopen("part.tmp", "w");
	  FPrintPartition(outf, best_part, 1);
	  fclose(outf);
	case 's':
	  outf = fopen("part.tmp", "w");
	  FPrintPartition(outf, best_part, 1);
	  fclose(outf);
	default:
	  break;
	}

	/* Update the temperature */
	T = T * Ts;

  } /* End of simulated annealing */

  // Complete mapping of the best partition on the network
  // (The previous ones where not updating the partition attributes)
  MapPartToNet(part, binet->net1);

  /* Free memory */
  free_d_mat(swwmat, nnod);
  free_d_vec(strength);
  RemovePartition(best_part);
  FreeLabelDict(nodeDict);
  free(glist);
  free(nlist);

  /* Done */
  return CompressPart(part);
}
