/*
  bipartite.c
  $LastChangedDate$
  $Revision$
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "graph.h"
#include "modules.h"

#include "bipartite.h"

#define EPSILON_MOD_B 1.e-6

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Network creation and removal
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Create an empty binet, whose networks point to NULL
  ---------------------------------------------------------------------
*/
struct binet *
CreateBipart()
{
  struct binet *temp;
  
  temp = (struct binet *)calloc(1, sizeof(struct binet));
  temp->net1 = NULL;
  temp->net2 = NULL;
  return temp;
}

/*
  ---------------------------------------------------------------------
  Create a bipartite network from a file.
  ---------------------------------------------------------------------
*/
struct binet *
FBuildNetworkBipart(FILE *inFile,
		    int weight_sw,
		    int add_weight_sw)
{
  struct node_gra *root1=NULL, *root2=NULL;
  void *dict1=NULL, *dict2=NULL;
  char label1[MAX_LABEL_LENGTH], label2[MAX_LABEL_LENGTH];
  struct node_gra *last1=NULL, *last2=NULL;
  struct node_gra *n1=NULL, *n2=NULL;
  struct binet *net=NULL;
  double weight;
  struct node_tree *nodeTree=NULL, *ntree1=NULL, *ntree2=NULL;
  int noReadItems = 0;
  char *line = NULL;
  size_t bufsiz = 0;
  ssize_t nbytes;

  /* Initialize the subnetworks */
  last1 = root1 = CreateHeaderGraph();
  last2 = root2 = CreateHeaderGraph();

  
  while ((nbytes = getline(&line, &bufsiz, inFile)) != -1){
    /* Read the labels (and weight, if necessary) */
    if (weight_sw == 0) {
      noReadItems = sscanf(line, "%s %s", label1, label2);
	  if(noReadItems != 2){
		printf ("Failed to read input: not enough fields in line %s (%d!=2). \n",line, noReadItems);
		return NULL;}
	  weight = 1;
    }
    else {
      noReadItems = sscanf(line,"%s %s %lf", label1, label2, &weight);
	  if(noReadItems != 3){
		printf ("Failed to read input: not enough fields in line %s (%d!=3). \n",line, noReadItems);
		return NULL;
	  }
    }

    /* Check if the nodes already exist, and create them otherwise */
    nodeTree = CreateNodeTree();
    strcpy(nodeTree->label, label1);
    ntree1 = *(struct node_tree **)tsearch((void *)nodeTree,
					   &dict1,
					   NodeTreeLabelCompare);
    if (ntree1->ref == NULL)
      ntree1->ref = last1 = CreateNodeGraph(last1, label1);
    else
      FreeNodeTree(nodeTree);
    n1 = ntree1->ref;

    nodeTree = CreateNodeTree();
    strcpy(nodeTree->label, label2);
    ntree2 = *(struct node_tree **)tsearch((void *)nodeTree,
					   &dict2,
					   NodeTreeLabelCompare);
    if (ntree2->ref == NULL)
      ntree2->ref = last2 = CreateNodeGraph(last2, label2);
    else
      FreeNodeTree(nodeTree);
    n2 = ntree2->ref;

    /* Create the links */
    AddAdjacency(n1, n2, 0, add_weight_sw, weight, 0);
    AddAdjacency(n2, n1, 0, add_weight_sw, weight, 0);
  }

  /* Free memory */
  FreeLabelDict(dict1);
  FreeLabelDict(dict2);
  free(line);

  /* Create the bipartite network and return */
  net = CreateBipart();
  net->net1 = root1;
  net->net2 = root2;
  return net;
}

/*
  ---------------------------------------------------------------------
  Build a modular bipartine network.

  modSizes: vector with the sizes of the modules (if NULL, then the
  modules are all of the same size nodpermod).

  nodpermod: number of nodes per module (NOTE: only used if
  modSizes==NULL)

  nmod: number of modules

  S2: number of 'teams'

  mmin (mmax): minimum (maximum) number of 'agents' per 'team'. Team
  size is determined by drawing a randomly generated integer uniformly
  distributed between mmin and mmax.

  p: probability of selecting agents of a given color for a team
  ---------------------------------------------------------------------
*/
struct binet *
BuildModularBipartiteNetwork(int *modSizes,
			     int nodpermod,
			     int nmod,
			     double *col_prob,
			     int S2,
			     int mmin, int mmax,
			     double geom_p,
			     double p,
			     gsl_rng *gen)
{
  int i, j;
  struct binet *net = NULL;
  struct node_gra *root1, *root2;
  struct node_gra **list1, **list2;
  struct node_gra *last_add1 = NULL, *last_add2 = NULL;
  int teamcol, target;
  int S1 = 0;
  int *module, *mod_starting, n_team_col;
  int m;
  int free_modSizes = 0;
  double dice, cum;
  char label[MAX_LABEL_LENGTH];
  
  /* Define modSizes if not defined yet */
  if (modSizes == NULL) {
    free_modSizes = 1; /* Remember to free memory at the end */
    modSizes = allocate_i_vec(nmod);
    for (i=0; i<nmod; i++) {
      modSizes[i] = nodpermod;
    }
  }

  /* Calculate the number of nodes S1 and allocate memory for the
     lists */
  S1 = 0;
  for (i=0; i<nmod; i++) {
    S1 += modSizes[i];
  }
  module = allocate_i_vec(S1);
  list1 = (struct node_gra **)calloc(S1, sizeof(struct node_gra *));
  list2 = (struct node_gra **)calloc(S2, sizeof(struct node_gra *));

  /* Assign each node to a module */
  S1 = 0;
  mod_starting = allocate_i_vec(nmod);
  for (i=0; i<nmod; i++) {
    mod_starting[i] = S1;  /* Number of the first node in the module */
    for (j=0; j<modSizes[i]; j++) {
      module[S1++] = i;
    }
  }

  /* Create the two networks */
  last_add1 = root1 = CreateHeaderGraph();
  last_add2 = root2 = CreateHeaderGraph();

  for (i = 0; i<S1; i++) {
    sprintf(label, "A%d", i + 1);
    last_add1 = list1[i] = CreateNodeGraph(last_add1, label);
    list1[i]->inGroup = module[i];
    list1[i]->ivar1 = 0;
  }
  for (i = 0; i<S2; i++) {
    sprintf(label, "B%d", i + 1);
    last_add2 = list2[i] = CreateNodeGraph(last_add2, label);
    list2[i]->ivar1 = 1;
  }

  /* Create the links */
  for (i = 0; i<S2; i++) {  /* loop over teams */

    /* Determine randomly which is the "color" of the team (according
       to the probabilities in col_prob if available, or with equal
       probability for all colors otherwise), that is, the module that
       will contribute, in principle, more nodes. */
    if (col_prob != NULL) {
      dice = gsl_rng_uniform(gen);
      cum = 0.0;
      teamcol = -1;
      while (cum < dice) {
	teamcol++;
	cum += col_prob[teamcol];
      }
    }
    else {
      teamcol = floor(gsl_rng_uniform(gen) * nmod);
    }
    list2[i]->inGroup = teamcol;
    n_team_col = 0; /* No nodes of the team color added so far */

    /* Determine the number of actors in the team */
    if (mmin < 0) { /* use a geometric distribution */
      m = geometric_dist_val(geom_p, gen);
    }
    else if (mmax != mmin) /* use uniform distribution */
      m = floor(gsl_rng_uniform(gen) * (double)(mmax+1-mmin) + mmin);
    else  /* all teams have the same size */
      m = mmin;

    for (j=0; j<m; j++) {  /* loop over spots in a team */

      if (gsl_rng_uniform(gen) < p &&
	  n_team_col < modSizes[teamcol]) { /* select node with the
					       module's color */
	n_team_col++;
	do {
	  target = mod_starting[teamcol] + floor(gsl_rng_uniform(gen) *
						 modSizes[teamcol]);
	} while (IsThereLink(list1[target], list2[i]) == 1);
      }
      else {  /* choose node at random */
	do {
	  target = floor(gsl_rng_uniform(gen) * S1);
	} while (IsThereLink(list1[target], list2[i]) == 1);
	if (module[target] == teamcol) {
	  n_team_col++;
	}
      }
      // Make the link
      AddAdjacency(list1[target], list2[i], 0, 0, 0, 0);
      AddAdjacency(list2[i], list1[target], 0, 0, 0, 0);
    } /* end of loop over spots */
  } /* end of loop over teams */
  
  /* Create the bipartite network */
  net = CreateBipart();
  net->net1 = root1;
  net->net2 = root2;

  /* Free memory */
  free(list1);
  free(list2);
  free_i_vec(module);
  free_i_vec(mod_starting);
  if (free_modSizes == 1) {
    free_i_vec(modSizes);
    modSizes = NULL;
  }

  /* Done */
  return net;
}

/*
  ---------------------------------------------------------------------
  Free the memory allocated to a bipartite network
  ---------------------------------------------------------------------
*/
void
RemoveBipart(struct binet *net)
{
  RemoveGraph(net->net1);
  RemoveGraph(net->net2);
  free(net);
  return;
}


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Node and multinode operations
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
/*
  ---------------------------------------------------------------------
  Calculates the number of common links between two nodes. The two
  nodes have to be in the same subnetwork; otherwise the subroutine
  returns 0.
  ---------------------------------------------------------------------
*/
int
NCommonLinksBipart(struct node_gra *n1, struct node_gra *n2)
{
  int ncom = 0;
  struct node_lis *l1 = n1->neig;
  struct node_lis *l2;

  while ((l1 = l1->next) != NULL) {
    l2 = n2->neig;
    while ((l2 = l2->next) != NULL) {
      if (l1->ref == l2->ref) ncom++;
    }
  }

  return ncom;
}

/*
  ---------------------------------------------------------------------
  Calculates Sum_a(w_ia*w_ja) for two nodes i and j. a is a node of the
  opposite subnetwork. The sum is over all nodes of that opposite subnetwork.
  The two nodes have to be in the same subnetwork; otherwise the subroutine
  returns 0.
  ---------------------------------------------------------------------
*/
double
SumProductsOfCommonWeightsBipart(struct node_gra *n1, struct node_gra *n2)
{
  double w_ia = 0.0, w_ja = 0.0, sww = 0.0;
  struct node_lis *l1 = n1->neig;
  struct node_lis *l2;

  while ((l1 = l1->next) != NULL) {
    w_ia = l1->weight;
    l2 = n2->neig;
    while ((l2 = l2->next) != NULL) {
      if (l1->ref == l2->ref) {
      	w_ja = l2->weight;
        sww += w_ia * w_ja;
      }
    }
  }
	
  return sww;
}

/*
  ---------------------------------------------------------------------
  Removes a node from the bipartite network. Variable set (1 or 2)
  indicates if the node lives in binet->net1 or binet->net2.
  ---------------------------------------------------------------------
*/
void
RemoveNodeBipart(struct binet *binet, char *label, int set)
{
  struct node_gra *p;
  struct node_gra *node;
  struct node_lis *nei;
  int count=0;

  /* Get to the target node */
  if (set == 1)
    p = binet->net1;
  else
    p = binet->net2;
  while (strcmp(p->next->label, label) != 0)
    p = p->next;
  node = p->next;

  /* Remove all links from/to the node */
  nei = node->neig;
  while (nei->next != NULL) {
    RemoveLink(node, nei->next->ref, 1);
    fprintf(stderr, "degree: %d\n", NodeDegree(node));
  }

  /* Remove the node */
  p->next = node->next;
  FreeNode(node);

  /* Renumber the nodes */
  if (set == 1)
    p = binet->net1;
  else
    p = binet->net2;
  while ((p = p->next) !=  NULL) {
    p->num = count++;
  }

  /* Renumber the node_lis nodes (adjacencies) */
  if (set == 1)
    p = binet->net2;
  else
    p = binet->net1;
  while ((p = p->next) !=  NULL) {
    nei = p->neig;
    while ((nei = nei->next) !=  NULL) {
      nei->node = nei->ref->num;
    }
  }

  return;
}

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Network operations
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Create a copy of a bipartite network
  ---------------------------------------------------------------------
*/
struct binet *
CopyBipart(struct binet *binet)
{
  struct binet *copy = NULL;
  struct node_gra *p1 = NULL, *p2 = NULL;
  struct node_gra *last1 = NULL, *last2 = NULL;
  struct node_lis *l = NULL;
  void *dict1 = NULL, *dict2 = NULL;
  struct node_tree *ntree1=NULL, *ntree2=NULL;

  /* Create the copy binet */
  copy = CreateBipart();
  last1 = copy->net1 = CreateHeaderGraph();
  last2 = copy->net2 = CreateHeaderGraph();

  /* Copy the nodes in each of the subnetworks and, simultaneously,
     create the label dictionaries for faster access to nodes */
  p1 = binet->net1;
  while ((p1 = p1->next) != NULL) {
    ntree1 = CreateNodeTree();
    strcpy(ntree1->label, p1->label);
    ntree1 = *(struct node_tree **)tsearch((void *)ntree1,
					   &dict1,
					   NodeTreeLabelCompare);
    ntree1->ref = last1 = CreateNodeGraph(last1, p1->label);
    last1->num = p1->num;
    last1->coorX = p1->coorX;
    last1->coorY = p1->coorY;
    last1->coorZ = p1->coorZ;
    last1->state = p1->state;
    last1->inGroup = p1->inGroup;
    last1->ivar1 = p1->ivar1;
    last1->dvar1 = p1->dvar1;
  }

  p2 = binet->net2;
  while ((p2 = p2->next) != NULL) {
    ntree2 = CreateNodeTree();
    strcpy(ntree2->label, p2->label);
    ntree2 = *(struct node_tree **)tsearch((void *)ntree2,
					   &dict2,
					   NodeTreeLabelCompare);
    ntree2->ref = last2 = CreateNodeGraph(last2, p2->label);
    last2->num = p2->num;
    last2->coorX = p2->coorX;
    last2->coorY = p2->coorY;
    last2->coorZ = p2->coorZ;
    last2->state = p2->state;
    last2->inGroup = p2->inGroup;
    last2->ivar1 = p2->ivar1;
    last2->dvar1 = p2->dvar1;
  }

  /* Copy the links */
  p1 = binet->net1;
  while ((p1 = p1->next) != NULL) {
    l = p1->neig;
    while ((l = l->next) != NULL) {
      AddAdjacency(GetNodeDict(p1->label, dict1),
		   GetNodeDict(l->nodeLabel, dict2),
		   0, 0, l->weight, 0);
      AddAdjacency(GetNodeDict(l->nodeLabel, dict2),
		   GetNodeDict(p1->label, dict1),
		   0, 0, l->weight, 0);
    }
  }

  /* Free memory */
  FreeLabelDict(dict1);
  FreeLabelDict(dict2);
  
  /* Done */
  return copy;
}

/*
  ---------------------------------------------------------------------
  Invert the networks in binet
  ---------------------------------------------------------------------
*/
struct binet *
InvertBipart(struct binet *net)
{
  struct node_gra *temp;

  temp = net->net1;
  net->net1 = net->net2;
  net->net2 = temp;

  return net;
}

/*
  ---------------------------------------------------------------------
  Counts the links in a bipartite network
  ---------------------------------------------------------------------
*/
int
NLinksBipart(struct binet *binet)
{
  struct node_gra *p = binet->net1;
  int nlink = 0;

  while ((p = p->next) !=  NULL)
    nlink +=  NodeDegree(p);

  return nlink;
}

/*
  ---------------------------------------------------------------------
  Project a binet into a weighted one mode network. The network is
  projected into the space of net1, and the weights in the projection
  are the number of common links in the binet.
  ---------------------------------------------------------------------
*/
struct node_gra *
ProjectBipart(struct binet *binet)
{
  struct node_gra *projnet = NULL;
  struct node_gra *last = NULL;
  struct node_gra *p1 = NULL, *p2 = NULL;
  void *dict = NULL;
  struct node_tree *ntree=NULL;
  int weight;

  /* Create the header of the projection network */
  last = projnet = CreateHeaderGraph();

  /* Create the nodes */
  p1 = binet->net1;
  while ((p1 = p1->next) != NULL) {
    ntree = CreateNodeTree();
    strcpy(ntree->label, p1->label);
    ntree = *(struct node_tree **)tsearch((void *)ntree,
					  &dict,
					  NodeTreeLabelCompare);
    ntree->ref = last = CreateNodeGraph(last, p1->label);
  }

  /* Create the links */
  p1 = binet->net1;
  while ((p1 = p1->next) != NULL) {
    p2 = p1;
    while ((p2 = p2->next) != NULL) {
      weight = NCommonLinksBipart(p1, p2);
      if (weight > 0) {
	AddAdjacency(GetNodeDict(p1->label, dict),
		     GetNodeDict(p2->label, dict),
		     0, 0, (double)weight, 0);
	AddAdjacency(GetNodeDict(p2->label, dict),
		     GetNodeDict(p1->label, dict),
		     0, 0, (double)weight, 0);
      }
    }
  }

  /* Free memory and return */
  FreeLabelDict(dict);
  return projnet;
}


struct node_gra *
ProjectBipartWeighted(struct binet *binet)
{
  struct node_gra *projnet = NULL;
  struct node_gra *last = NULL;
  struct node_gra *p1 = NULL, *p2 = NULL;
  void *dict = NULL;
  struct node_tree *ntree=NULL;
  int weight;

  /* Create the header of the projection network */
  last = projnet = CreateHeaderGraph();

  /* Create the nodes */
  p1 = binet->net1;
  while ((p1 = p1->next) != NULL) {
    ntree = CreateNodeTree();
    strcpy(ntree->label, p1->label);
    ntree = *(struct node_tree **)tsearch((void *)ntree,
					  &dict,
					  NodeTreeLabelCompare);
    ntree->ref = last = CreateNodeGraph(last, p1->label);
  }

  /* Create the links */
  p1 = binet->net1;
  while ((p1 = p1->next) != NULL) {
    p2 = p1;
    while ((p2 = p2->next) != NULL) {
	  weight = (double)SumProductsOfCommonWeightsBipart(p1,p2);
      if (weight > 0) {
	AddAdjacency(GetNodeDict(p1->label, dict),
		     GetNodeDict(p2->label, dict),
		     0, 0, (double)weight, 0);
	AddAdjacency(GetNodeDict(p2->label, dict),
		     GetNodeDict(p1->label, dict),
		     0, 0, (double)weight, 0);
      }
    }
  }

  /* Free memory and return */
  FreeLabelDict(dict);
  return projnet;
}


/*
  ---------------------------------------------------------------------
  Randomize links in a bipartite network
  ---------------------------------------------------------------------
*/
struct binet *
RandomizeBipart(struct binet *binet, double times, gsl_rng *gen)
{
  int i;
  int nlink, niter, coun = 0;
  int target1, target2;
  struct node_gra *p = NULL;
  struct node_lis *l = NULL;
  struct node_gra *n1, *n2, *n3, *n4;
  struct node_gra **ori, **des;

  /* Build the link lists (one for link origins and one for ends) */
  nlink = NLinksBipart(binet);
  niter = ceil(times * (double)nlink);
  ori = (struct node_gra **) calloc(nlink, sizeof(struct node_gra *));
  des = (struct node_gra **) calloc(nlink, sizeof(struct node_gra *));

  p = binet->net1;
  while ((p = p->next) !=  NULL) {
    l = p->neig;
    while ((l = l->next) !=  NULL) {
      ori[coun] = p;
      des[coun] = l->ref;
      coun++;
    }
  }

  if (coun !=  nlink)
    fprintf(stderr, "Error in RandomizeBipart: coun !=  nlink!!\n");

  /* Randomize the network */
  for (i=0; i<niter; i++) {
    
    /* select the two links (four nodes) to swap */
    do {
      target1 = floor(gsl_rng_uniform(gen) * (double)nlink);
      n1 = ori[target1];
      n2 = des[target1];

      do {
	target2 = floor(gsl_rng_uniform(gen) * (double)nlink);
	n3 = ori[target2];
	n4 = des[target2];
      } while (n1 == n3 || n2 == n4);

    } while (IsThereLink(n1, n4) == 1 ||
	     IsThereLink(n2, n3) == 1);

    /* switch the link */
    RemoveLink(n1, n2, 1);
    RemoveLink(n3, n4, 1);
    AddAdjacency(n1, n4, 0, 0, 0, 1);
    AddAdjacency(n4, n1, 0, 0, 0, 1);
    AddAdjacency(n3, n2, 0, 0, 0, 1);
    AddAdjacency(n2, n3, 0, 0, 0, 1);

    ori[target1] = n1;
    des[target1] = n4;
    ori[target2] = n3;
    des[target2] = n2;
  }

  free(ori);
  free(des);

  return binet;
}


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Network output
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Print a bipartite network in Pajek format. Node coordinates are
  printed if coor_sw is not 0. Weight is printed if weight_sw is not
  0. If symmetric_sw == 1, the link is only printed once (a-b),
  otherwise it will be printed twice if necessary (a-b, and b-a).
  ---------------------------------------------------------------------
*/
void
FPrintPajekFileBipart (char *fname,
		       struct binet *binet,
		       int coor_sw,
		       int weight_sw)
{
  struct node_gra *p=NULL;
  struct node_lis *n=NULL;
  FILE *outF;
  int Stot, S1, S2;

  /* Preliminaries */
  S1 = CountNodes(binet->net1);
  S2 = CountNodes(binet->net2);
  Stot  = S1 + S2;

  outF = fopen(fname, "w");

  fprintf(outF, "*Vertices %d %d\n", Stot, S1);

  /* Print the nodes */
  p = binet->net1;
  while((p = p->next) != NULL) {
    if (coor_sw == 1) {
      fprintf(outF, "%d \"%s\"           %lf %lf %lf\n",
	      p->num + 1, p->label, p->coorX, p->coorY, p->coorZ);
    }
    else {
      fprintf(outF, "%d \"%s\"\n", p->num + 1, p->label);
    }
  }

  p = binet->net2;
  while((p = p->next) != NULL) {
    if (coor_sw == 1) {
      fprintf(outF, "%d \"%s\"           %lf %lf %lf\n",
	      p->num + S1 + 1, p->label, p->coorX, p->coorY, p->coorZ);
    }
    else {
      fprintf(outF, "%d \"%s\"\n", p->num + S1 + 1, p->label);
    }
  }

  /* Print the links title */
  fprintf(outF,"*Edges\n");

  /* Print the links */
  p = binet->net1;
  while ((p = p->next) != NULL) {
    n = p->neig;
    while ((n = n->next) != NULL) {
      if (weight_sw == 0)
	fprintf(outF, "%d   %d\n",
		p->num + 1, n->ref->num + S1 + 1);
      else
	fprintf(outF, "%d   %d   %g\n", 
		p->num + 1, n->ref->num + S1 + 1, n->weight);
    }
  }

  /* Done */
  fclose(outF);
}

/*
  ---------------------------------------------------------------------
  Print a bipartite network in two-column format.
  ---------------------------------------------------------------------
*/
void
FPrintBipart (FILE *outf, struct binet *binet, int weight_sw)
{
  struct node_gra *p=NULL;
  struct node_lis *n=NULL;

  /* Print the links */
  p = binet->net1;
  while ((p = p->next) != NULL) {
    n = p->neig;
    while ((n = n->next) != NULL) {
      if (weight_sw == 0)
        fprintf(outf,
                "%s %s\n", 
                p->label, n->ref->label);
      else
        fprintf(outf,
                "%s %s %g\n", 
                p->label, n->ref->label, n->weight);

    }
  }
}

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Network modularity
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Calculate the modularity of a partition of a bipartite network
  ---------------------------------------------------------------------
*/
double
ModularityBipart(struct binet *binet, struct group *part)
{
  struct node_gra *p = binet->net2;
  struct node_lis *n1, *n2;
  double bimod = 0.0;
  int c12;
  int t1, t2;
  double sms = 0.0, s2ms = 0.0, sms2 = 0.0;

  /* Calculate sms=sum(m_s) and sms2=sum(m_s^2) */
  while ((p = p->next) != NULL) {
    sms += (double)NodeDegree(p);
    sms2 += (double)(NodeDegree(p) * NodeDegree(p));
  }
  s2ms = sms * sms;

  /* Calculate the modularity */
  while ((part = part->next) != NULL){

    n1 = part->nodeList;
    while ((n2 = n1 = n1->next) != NULL) {
      while ((n2 = n2->next) != NULL) {
        c12 = NCommonLinksBipart(n1->ref, n2->ref);
        t1 = NodeDegree(n1->ref);
        t2 = NodeDegree(n2->ref);

	      bimod += (double)c12 / (sms2 - sms) - (double)(t1 * t2) / s2ms;
      }
    }
  }
  bimod *= 2.;

  /* Done */
  return bimod;
}

/*
  ---------------------------------------------------------------------
  Calculate the modularity of a partition of a weighted bipartite network
  ---------------------------------------------------------------------
*/
double
ModularityBipartWeighted(struct binet *binet, struct group *part)
{
  struct node_gra *p = binet->net2;
  struct node_lis *n1, *n2;
  double bimod = 0.0;
  double sww12;
  double s1, s2;
  double sWa = 0.0, s2Wa = 0.0, sWa2 = 0.0;

  /* Calculate s2Wa=(sum W_a)^2 and sWa2=sum(W_a^2) */
  while ((p = p->next) != NULL) {
    sWa += (double)NodeStrength(p);
    sWa2 += (double)(NodeStrength(p) * NodeStrength(p));
  }
  s2Wa = sWa * sWa;

  /* Calculate the modularity */
  while ((part = part->next) != NULL){

    n1 = part->nodeList;
    while ((n2 = n1 = n1->next) != NULL) {
      while ((n2 = n2->next) != NULL) {
        sww12 = (double)SumProductsOfCommonWeightsBipart(n1->ref, n2->ref);
        s1 = (double)NodeStrength(n1->ref);
        s2 = (double)NodeStrength(n2->ref);
        
        bimod += sww12 / (sWa2 - sWa) - (s1 * s2) / s2Wa;
      }
    }
  }
  bimod *= 2.;

  /* Done */
  return bimod;
}

/*
  ---------------------------------------------------------------------
  Calculate the modularity of a partition of a weighted bipartite network
  ---------------------------------------------------------------------
*/
double ModularityBipartWeightedFast(struct binet *binet, struct group *part, double **swwmat)
{
  struct node_lis *n1, *n2;
  double bimod = 0.0;

  /* Calculate the modularity */
  while ((part = part->next) != NULL){
    n1 = part->nodeList;
    while ((n2 = n1 = n1->next) != NULL) {
      while ((n2 = n2->next) != NULL) {
        bimod += swwmat[n1->node][n2->node];
      }
    }
  }

  /* Done */
  return bimod;
}



/**
Given a bipartite network and a module partition, compute the
bipartite role on the projection of the first component and
write the output in outf using the following format:

Label \t Module \t Role \t P \ z

Example: 
Mynode\t1\tR3\t0.6500\t-1.4400
**/
void
FPrintTabNodesBipart(FILE *outf, struct binet *network,  struct group *modules, int degree_based,int weighted)
{
  struct node_gra *projected;
  struct group    *g = NULL;
  struct node_lis *n = NULL;
  double P, z;
  int group_nb = 0;
  int role;

  // Project the first component of the bipartite network.
  if (weighted == 0 || degree_based == 1)
	projected = ProjectBipart(network);
  else
	projected = ProjectBipartWeighted(network);
  
  // Loop through the modules of the bipartite network and compute the
  // roles on the projection.
  MapPartToNetFast(modules, projected);
  g = modules;
  while ((g = g->next) != NULL) {
    n = g->nodeList;
	group_nb ++;
    while ((n = n->next) != NULL) {
	  if (degree_based==1){
		P = ParticipationCoefficient(n->ref);
		z = WithinModuleRelativeDegree(n->ref, g);
		role = GetRole(P,z) + 1;
	  }
	  else{
		P = WeightedParticipationCoefficient(n->ref,modules);
		z = WithinModuleRelativeStrength(n->ref, g);
		role = GetRole(P,z) + 1;
	  }
	  printf ("%-20s\t%d\tR%d\t%f\t%f\n",
			  n->nodeLabel,
			  group_nb,
			  role,
			  P,
			  z);
	}}
  // Remap the modules on the original network. 
  MapPartToNet(modules, network->net1);

  // Free Memory 
  RemoveGraph(projected);
}



/*
  ---------------------------------------------------------------------
  Calculte the participation coefficient of a node in a bipartite
  network
  ---------------------------------------------------------------------
*/
double
ParticipationCoefficientBipart(struct node_gra *node)
{
  struct node_lis *aTeamPoint = node->neig, *aNodePoint = NULL;
  struct node_gra *temp = node;
  int norm = 0;
  int maxNum = 0;
  int *inGroupCount = NULL, i;
  double P = 0.0;

  /* Determine the maximum number of groups and create an array for
     inGroup counts */
  while (temp != NULL) {
    maxNum = temp->num;
    temp = temp->next;
  }
  inGroupCount = allocate_i_vec(maxNum + 1);
  for (i=0; i<=maxNum; i++)
    inGroupCount[i] = 0;

  /* Go through the teams to which the node belongs */
  while ((aTeamPoint = aTeamPoint->next) != NULL) {
    /* Go through the members of the team and update inGroupCount */
    aNodePoint = aTeamPoint->ref->neig;
    while ((aNodePoint = aNodePoint->next) != NULL) {
      if (aNodePoint->ref != node) {
	norm++;
	inGroupCount[aNodePoint->ref->inGroup] += 1;
      }
    }
  }

  /* Compute the participation coefficient */
  for (i=0; i<=maxNum; i++)
    P +=
      (double)inGroupCount[i] * (double)inGroupCount[i] /
      ((double)norm * (double)norm);
  P = 1.0 - P;

  /* Free memory */
  free_i_vec(inGroupCount);

  /* Done */
  return P;
}

/*
  ---------------------------------------------------------------------
  Calculate the largest, smallest, average, and standard deviation of
  the participation coefficient of all nodes in a network.
  ---------------------------------------------------------------------
*/
void
StatisticsParticipationCoefficientBipart(struct node_gra *net,
					 double *PMean,
					 double *PStddev,
					 double *PMin,
					 double *PMax
					 )
{
  struct node_gra *p = net;
  int N = CountNodes(net);
  double *PList = NULL;
  
  /* Allocate memory */
  PList = allocate_d_vec(N);

  /* Calculate all Ps */
  while ((p = p->next) != NULL)
    PList[p->num] = ParticipationCoefficientBipart(p);
  
  /* Get the statistics */
  *PMean = mean(PList, N);
  *PStddev = stddev(PList, N);
  *PMin = min(PList, N);
  *PMax = max(PList, N);

  /* Free memory */
  free_d_vec(PList);

  /* Done */
  return;
}

/* /\* */
/*   --------------------------------------------------------------------- */
/*   Calculate the modularity of a co-partition of both subnetworks of a */
/*   bipartite network */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* double BipartCoModularity(struct binet *binet, struct group *part) */
/* { */
/*   struct node_lis *n1, *n2; */
/*   double bimod = 0.0; */
/*   int t1, t2; */
/*   int L; */
/* /\*   int S1, S2; *\/ */

/* /\*   S1 = CountNodes(binet->net1); *\/ */
/* /\*   S2 = CountNodes(binet->net2); *\/ */
/*   L = NLinksBipart(binet); */

/*   // Calculate the modularity */
/*   while ((part = part->next) != NULL){ */
/*     n1 = part->nodeList; */
/*     while ((n2 = n1 = n1->next) != NULL) { */
/*       while ((n2 = n2->next) != NULL) { */
/* 	if ( n1->ref->ivar1 != n2->ref->ivar1 ) { */
/* 	  t1 = NodeDegree(n1->ref); */
/* 	  t2 = NodeDegree(n2->ref); */
/* 	  bimod += (IsThereLink(n1->ref, n2->node) -  */
/* 		    (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	} */
/*       } */
/*     } */
/*   } */
  
/*   return bimod; */
/* } */





/* /\* */
/*   --------------------------------------------------------------------- */
/*   Build a bipartite network with the nodes in a certain module */
/*   only. The module can contain nodes in both subnetworks of the */
/*   bipartite network. */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* struct binet *BuildNetworkFromCopart(struct group *copart) */
/* { */
/*   struct binet *binet = NULL; */
/*   struct node_gra *root1 = NULL; */
/*   struct node_gra *root2 = NULL; */
/*   struct node_gra *list1[maxim_int], *ref1[maxim_int]; */
/*   struct node_gra *list2[maxim_int], *ref2[maxim_int]; */
/*   struct node_gra *last_add1 = NULL; */
/*   struct node_gra *last_add2 = NULL; */
/*   struct node_lis *node = NULL; */
/*   struct node_gra *n1 = NULL, *n2 = NULL; */
/*   int S1 = 0, S2 = 0;  */
/*   int i, j; */

/*   // Create the subnetworks */
/*   last_add1 = root1 = CreateHeaderGraph(); */
/*   last_add2 = root2 = CreateHeaderGraph(); */

/*   // Add the nodes */
/*   node = copart->nodeList; */
/*   while ((node = node->next) != NULL) { */
/*     if (node->ref->ivar1 == 0) { */
/*       ref1[S1] = node->ref; */
/*       last_add1 = list1[S1] = CreateNodeGraph(last_add1,  */
/* 					      node->node, */
/* 					      node->node, */
/* 					      node->node); */
/*       list1[S1]->ivar1 = 0; */
/*       S1++; */
/*     } */
/*     else { */
/*       ref2[S2] = node->ref; */
/*       last_add2 = list2[S2] = CreateNodeGraph(last_add2,  */
/* 					      node->node, */
/* 					      node->node, */
/* 					      node->node); */
/*       list2[S2]->ivar1 = 1; */
/*       S2++; */
/*     } */
/*   } */

/*   // Add the links */
/*   for (i=0; i<S1; i++) { */
/*     for (j=0; j<S2; j++) { */
/*       if (IsThereLink(ref1[i], ref2[j]->num) == 1) { */
/* 	AddAdjacencyFull(list1[i], list2[j], 10); */
/* 	AddAdjacencyFull(list2[j], list1[i], 10); */
/* 	GetLink(list1[i], list2[j]->num)->weight =  */
/* 	  GetLink(ref1[i], ref2[j]->num)->weight; */
/* 	GetLink(list2[j], list1[i]->num)->weight =  */
/* 	  GetLink(ref2[j], ref1[i]->num)->weight; */
/*       } */
/*     } */
/*   } */
  
/*   // Create the bipartite network and return it */
/*   binet = CreateBipart(); */
/*   binet->net1 = root1; */
/*   binet->net2 = root2; */
/*   return binet; */
/* } */


/* /\* */
/*   --------------------------------------------------------------------- */
/*   Splits actors and teams in one group into two new groups */
/*   (thermalized with the outer SA) */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* ThermalBipartworkCoSplit(struct group *target_g, struct group *empty_g, */
/* 			double Ti, double Tf, */
/* 			gsl_rng *gen) */
/* { */
/*   struct group *glist[2]; */
/*   struct node_gra *nlist[maxim_int]; */
/*   struct node_gra *nlist_ori0[maxim_int]; */
/*   struct node_gra *nlist_ori1[maxim_int]; */
/*   struct node_gra *p = NULL; */
/*   int nnod = 0; */
/*   int i, j; */
/*   int target, oldg, newg; */
/*   double dE = 0.0; */
/*   double T, Ts = 0.95; */
/*   double dice; */
/*   struct binet *module_binet = NULL; */
/*   struct group *temp_part = NULL; */
/*   struct node_lis *nod = NULL; */
/*   int S1, S2, t1, t2, L; */
/*   int *nlinks = NULL; */
/*   int **isThereLink = NULL; */

/*   // Create a bipartite network with the nodes in the target_g only */
/*   module_binet = BuildNetworkFromCopart(target_g); */
/*   L = NLinksBipart(module_binet); */
/*   S1 = CountNodes(module_binet->net1); */
/*   S2 = CountNodes(module_binet->net2); */
/*   nlinks = allocate_i_vec(S1 + S2); */
/*   isThereLink = allocate_i_mat(S1 + S2, S1 + S2); */
/*   // Normalize the link matrix */
/*   for (i=0; i<S1+S2; i++) { */
/*     for (j=0; j<S1+S2; j++) { */
/*       isThereLink[i][j] = isThereLink[j][i] = 0; */
/*       /\* 	IsThereLink(nlist[i], nlist[j]->num); *\/ */
/*     } */
/*   } */

/*   // Create the groups */
/*   temp_part = CreateHeaderGroup(); */
/*   glist[0] = CreateGroup(temp_part, 0); */
/*   glist[1] = CreateGroup(temp_part, 1); */
  
/*   // Randomly assign the nodes to the groups and count their links */
/*   p = module_binet->net1; */
/*   while ((p = p->next) != NULL) { */
/*     p->trans = nnod++; */
/*     nlist[p->trans] = p; */
/*     nlinks[p->trans] = NodeDegree(p); */
/*     dice = gsl_rng_uniform(gen); */
/*     if (dice < 0.5) { */
/*       AddNodeToGroup(glist[0], p); */
/*     } */
/*     else { */
/*       AddNodeToGroup(glist[1], p); */
/*     } */
/*   } */
/*   p = module_binet->net2; */
/*   while ((p = p->next) != NULL) { */
/*     p->trans = nnod++; */
/*     nlist[p->trans] = p; */
/*     nlinks[p->trans] = NodeDegree(p); */
/*     dice = gsl_rng_uniform(gen); */
/*     if (dice < 0.5) { */
/*       AddNodeToGroup(glist[0], p); */
/*     } */
/*     else { */
/*       AddNodeToGroup(glist[1], p); */
/*     } */
/*     // Find the links */
/*     nod = p->neig; */
/*     while ((nod = nod->next) != NULL) { */
/*       isThereLink[p->trans][nod->ref->trans] =  */
/* 	isThereLink[nod->ref->trans][p->trans] = 1; */
/*     } */
/*   } */

/*   // Do SA to "optimize" the splitting */
/*   T = Ti; */
/*   while (T >= Tf) { */
    
/*     for (i=0; i<nnod; i++) { */
/*       target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       newg = 1 - oldg; */
	
/*       // Calculate the change of energy */
/*       dE = 0.0; */
/*       t1 = nlinks[nlist[target]->trans]; */

/*       // Old group contribution */
/*       nod = glist[oldg]->nodeList; */
/*       while ((nod = nod->next) != NULL) { */
/* 	if (nlist[target]->ivar1 != nod->ref->ivar1) { */
/* 	  t2 = nlinks[nod->ref->trans]; */
/* 	  dE -= (isThereLink[target][nod->ref->trans] - */
/* 		 (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	} */
/*       } */

/*       // New group contribution */
/*       nod = glist[newg]->nodeList; */
/*       while ((nod = nod->next) != NULL) { */
/* 	if (nlist[target]->ivar1 != nod->ref->ivar1) { */
/* 	  t2 = nlinks[nod->ref->trans]; */
/* 	  dE += (isThereLink[target][nod->ref->trans] - */
/* 		 (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	} */
/*       } */
      
/*       // Accept the change according to the Boltzman factor */
/*       if( (dE >= 0.0) || (gsl_rng_uniform(gen) < exp(dE/T)) ){ */
/* 	MoveNode(nlist[target], glist[oldg], glist[newg]); */
/*       } */

/*     } // End of for loop */

/*     T = T * Ts; */
/*   } // End of SA */

/*   // Actually move the nodes in the original network */
/*   // List the nodes in the original network */
/*   nod = target_g->nodeList; */
/*   while ((nod = nod->next) != NULL) { */
/*     if (nod->ref->ivar1 == 0) */
/*       nlist_ori0[nod->node] = nod->ref; */
/*     else */
/*       nlist_ori1[nod->node] = nod->ref; */
/*   } */
/*   // Move nodes around */
/* /\*   printf("\t---\n"); *\/ */
/*   nod = glist[1]->nodeList; */
/*   while ((nod = nod->next) != NULL) { */
/*     if (nod->ref->ivar1 == 0) { */
/*       MoveNode(nlist_ori0[nod->node], target_g, empty_g); */
/*     } */
/*     else { */
/*       MoveNode(nlist_ori1[nod->node], target_g, empty_g); */
/*     } */
/*   } */

/*   // Free memory */
/*   RemoveBipart(module_binet); */
/*   RemovePartition(temp_part); */
/*   free_i_vec(nlinks); */
/*   free_i_mat(isThereLink, S1 + S2); */
/* } */



/* /\* */
/*   --------------------------------------------------------------------- */
/*   Create a partition according to the isolated clusters in the network */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* SetClusterInGroup(struct node_gra *node, int *ngroup) */
/* { */
/*   struct node_lis *nei; */

/*   node->inGroup = *ngroup; */
/*   nei = node->neig; */
/*   while ((nei = nei->next) != NULL) { */
/*     if (nei->ref->inGroup < 0) { */
/*       SetClusterInGroup(nei->ref, ngroup); */
/*     } */
/*   } */
/* } */
/* struct group *BipartClustersPartition(struct binet *binet) */
/* { */
/*   struct node_gra *p; */
/*   int ngroup = 0, i; */
/*   struct group *part = CreateHeaderGroup(); */
/*   struct group **glist, *lastg; */

/*   ResetNetGroup(binet->net1); */
/*   ResetNetGroup(binet->net2); */

/*   // Set the inGroup of each node according to the cluster it belongs */
/*   // to */
/*   p = binet->net1; */
/*   while ((p = p->next) != NULL) { */
/*     if (p->inGroup < 0) { */
/*       SetClusterInGroup(p, &ngroup); */
/*       ngroup++; */
/*     } */
/*   } */
/*   p = binet->net2; */
/*   while ((p = p->next) != NULL) { */
/*     if (p->inGroup < 0) { */
/*       SetClusterInGroup(p, &ngroup); */
/*       ngroup++; */
/*     } */
/*   } */

/*   // Create the groups */
/*   glist = (struct group **)malloc(ngroup * sizeof(struct group *)); */
/*   lastg = part; */
/*   for (i=0; i<ngroup; i++) {  */
/*     lastg = glist[i] = CreateGroup(lastg, i); */
/*   } */

/*   // Add the nodes to the groups */
/*   p = binet->net1; */
/*   while ((p = p->next) != NULL) { */
/*     AddNodeToGroup(glist[p->inGroup], p); */
/*   } */
/*   p = binet->net2; */
/*   while ((p = p->next) != NULL) { */
/*     AddNodeToGroup(glist[p->inGroup], p); */
/*   } */
  
/*   return part; */
/* } */


/* /\* */
/*   --------------------------------------------------------------------- */
/*   Same as ThermalBipartCoSplit, but checks first if there are */
/*   disconnected components. If indeed the module is disconnected, */
/*   proposes a split using this fact. */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* ThermalPercBipartworkCoSplit(struct group *target_g, struct group *empty_g, */
/* 			    double prob, double Ti, double Tf, */
/* 			    gsl_rng *gen) */
/* { */
/*   struct group *glist[maxim_int], *g = NULL; */
/*   struct node_gra *nlist[maxim_int]; */
/*   struct node_gra *nlist_ori0[maxim_int]; */
/*   struct node_gra *nlist_ori1[maxim_int]; */
/*   struct node_gra *p = NULL; */
/*   int nnod = 0; */
/*   int i, j; */
/*   int target, oldg, newg; */
/*   double dE = 0.0; */
/*   double T, Ts = 0.95; */
/*   double dice; */
/*   struct binet *module_binet = NULL; */
/*   struct group *temp_part = NULL; */
/*   struct node_lis *nod = NULL; */
/*   int S1, S2, t1, t2, L; */
/*   int *nlinks = NULL; */
/*   int **isThereLink = NULL; */
/*   int ngroup, nmerge; */
/*   int node1, mod1, node2, mod2; */

/*   // Create a bipartite network with the nodes in the target_g only */
/*   module_binet = BuildNetworkFromCopart(target_g); */
/*   L = NLinksBipart(module_binet); */
/*   S1 = CountNodes(module_binet->net1); */
/*   S2 = CountNodes(module_binet->net2); */

/*   // Map the nodes to a the node lists for faster access */
/*   nnod = 0; */
/*   p = module_binet->net1; */
/*   while ((p = p->next) != NULL) { */
/*     p->trans = nnod++; */
/*     nlist[p->trans] = p; */
/*   } */
/*   p = module_binet->net2; */
/*   while ((p = p->next) != NULL) { */
/*     p->trans = nnod++; */
/*     nlist[p->trans] = p; */
/*   } */

/*   // Check if there are disconnected clusters */
/*   temp_part = BipartClustersPartition(module_binet); */

/*   if (CountGroups(temp_part) > 1 && gsl_rng_uniform(gen) < prob) { */
/*     /\* */
/*       Determine modules using clusters */
/*     *\/ */
/*     // Map the groups to an array for faster access */
/*     ngroup = 0; */
/*     g = temp_part; */
/*     while ((g = g->next) != NULL) { */
/*       if (ngroup != g->label) printf("ERROR\n"); */
/*       glist[ngroup++] = g; */
/*     } */

/*     // Merge pairs of modules until only two are left */
/*     nmerge = ngroup-2; */
/*     for (i=0; i<nmerge; i++) { */
/*       node1 = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       mod1 = nlist[node1]->inGroup; */
/*       do { */
/* 	node2 = floor(gsl_rng_uniform(gen) * (double)nnod); */
/* 	mod2 = nlist[node2]->inGroup; */
/*       } while (mod2 == mod1); */
/*       MergeGroups(glist[mod1], glist[mod2]); */
/*     } */
/*     CompressPart(temp_part); */
/*     glist[0] = temp_part->next; */
/*     glist[1] = temp_part->next->next; */
/*   } */
/*   else { */
/*     /\* */
/*       Determine modules using SA */
/*     *\/ */
/*     // Create some matrices for faster access to the number of links */
/*     // and to the links themselves */
/*     nlinks = allocate_i_vec(S1 + S2); */
/*     isThereLink = allocate_i_mat(S1 + S2, S1 + S2); */
/*     for (i=0; i<S1+S2; i++) { */
/*       for (j=0; j<S1+S2; j++) { */
/* 	isThereLink[i][j] = isThereLink[j][i] = 0; */
/* 	/\* 	IsThereLink(nlist[i], nlist[j]->num); *\/ */
/*       } */
/*     } */
    
/*     // Create the groups */
/*     RemovePartition(temp_part); */
/*     temp_part = CreateHeaderGroup(); */
/*     glist[0] = CreateGroup(temp_part, 0); */
/*     glist[1] = CreateGroup(temp_part, 1); */
    
/*     // Randomly assign the nodes to the groups and count their links */
/*     p = module_binet->net1; */
/*     while ((p = p->next) != NULL) { */
/*       nlinks[p->trans] = NodeDegree(p); */
/*       dice = gsl_rng_uniform(gen); */
/*       if (dice < 0.5) { */
/* 	AddNodeToGroup(glist[0], p); */
/*       } */
/*       else { */
/* 	AddNodeToGroup(glist[1], p); */
/*       } */
/*     } */
/*     p = module_binet->net2; */
/*     while ((p = p->next) != NULL) { */
/*       nlinks[p->trans] = NodeDegree(p); */
/*       dice = gsl_rng_uniform(gen); */
/*       if (dice < 0.5) { */
/* 	AddNodeToGroup(glist[0], p); */
/*       } */
/*       else { */
/* 	AddNodeToGroup(glist[1], p); */
/*       } */
/*       // Find the links */
/*       nod = p->neig; */
/*       while ((nod = nod->next) != NULL) { */
/* 	isThereLink[p->trans][nod->ref->trans] =  */
/* 	  isThereLink[nod->ref->trans][p->trans] = 1; */
/*       } */
/*     } */
    
/*     // Do SA to "optimize" the splitting */
/*     T = Ti; */
/*     while (T >= Tf) { */
      
/*       for (i=0; i<nnod; i++) { */
/* 	target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/* 	oldg = nlist[target]->inGroup; */
/* 	newg = 1 - oldg; */
	
/* 	// Calculate the change of energy */
/* 	dE = 0.0; */
/* 	t1 = nlinks[nlist[target]->trans]; */
	
/* 	// Old group contribution */
/* 	nod = glist[oldg]->nodeList; */
/* 	while ((nod = nod->next) != NULL) { */
/* 	  if (nlist[target]->ivar1 != nod->ref->ivar1) { */
/* 	    t2 = nlinks[nod->ref->trans]; */
/* 	    dE -= (isThereLink[target][nod->ref->trans] - */
/* 		   (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	  } */
/* 	} */
	
/* 	// New group contribution */
/* 	nod = glist[newg]->nodeList; */
/* 	while ((nod = nod->next) != NULL) { */
/* 	  if (nlist[target]->ivar1 != nod->ref->ivar1) { */
/* 	    t2 = nlinks[nod->ref->trans]; */
/* 	    dE += (isThereLink[target][nod->ref->trans] - */
/* 		   (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	  } */
/* 	} */
	
/* 	// Accept the change according to the Boltzman factor */
/* 	if( (dE >= 0.0) || (gsl_rng_uniform(gen) < exp(dE/T)) ){ */
/* 	  MoveNode(nlist[target], glist[oldg], glist[newg]); */
/* 	} */
	
/*       } // End of for loop */
      
/*       T = T * Ts; */
/*     } // End of SA */

/*     // Free some memory */
/*     free_i_mat(isThereLink, S1 + S2); */
/*   } */

/*   /\* */
/*     Regardless of the procedure used, glist[0], glist[1] contain now */
/*     the suggested partition. */
/*   *\/ */
/*   // Actually move the nodes in the original network */
/*   // List the nodes in the original network */
/*   nod = target_g->nodeList; */
/*   while ((nod = nod->next) != NULL) { */
/*     if (nod->ref->ivar1 == 0) */
/*       nlist_ori0[nod->node] = nod->ref; */
/*     else */
/*       nlist_ori1[nod->node] = nod->ref; */
/*   } */
/*   // Move nodes around */
/*   nod = glist[1]->nodeList; */
/*   while ((nod = nod->next) != NULL) { */
/*     if (nod->ref->ivar1 == 0) { */
/*       MoveNode(nlist_ori0[nod->node], target_g, empty_g); */
/*     } */
/*     else { */
/*       MoveNode(nlist_ori1[nod->node], target_g, empty_g); */
/*     } */
/*   } */

/*   // Free memory */
/*   RemoveBipart(module_binet); */
/*   RemovePartition(temp_part); */
/*   free_i_vec(nlinks); */
/* } */




/* /\* */
/*   --------------------------------------------------------------------- */
/*   Identify modules simultaneously in both "subnetworks" of a bipartite */
/*   network using simulated annealing. */

/*   Ti: Initial temperature for the SA */
/*   Tf: Final temperature for the SA */
/*   Ts: Cooling factor */
/*   fac: Iteration factor */
/*   merge: implement collective moves (merge=1) or not (merge=0) */
/*   prob: only used if merge==1. Use percolation split with probability */
/*   prob (prob>0) or do not use percolation (prob<0) */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* struct group *SACommunityCoIdentBipart(struct binet *binet, */
/* 				      double Ti, double Tf, double Ts, */
/* 				      double fac, */
/* 				      int merge, */
/* 				      double prob, */
/* 				      gsl_rng *gen) */
/* { */
/*   int i; */
/*   struct node_gra *net1 = binet->net1; */
/*   struct node_gra *net2 = binet->net2; */
/*   struct node_gra *p, *p2; */
/*   int S1, S2, L, nnod; */
/*   struct group *part = NULL, *g = NULL, *split = NULL; */
/*   struct node_gra *nlist[maxim_int]; */
/*   struct group *glist[maxim_int], *lastg; */
/*   int cicle1, cicle2; */
/*   int count = 0, limit = 25; */
/*   double energy, energyant, dE, e; */
/*   double T; */
/*   int target, oldg, newg; */
/*   int g1, g2, empty; */
/*   struct node_lis *nod, *nod2; */
/*   double t1, t2; */
/*   struct group *best_part = NULL; */
/*   double best_E = -100.0; */

/*   // Count nodes and links */
/*   S1 = CountNodes(binet->net1); */
/*   S2 = CountNodes(binet->net2); */
/*   L = NLinksBipart(binet); */
/*   nnod = S1 + S2; */

/*   // Create the groups and assign each node to one group */
/*   ResetNetGroup(net1); // All nodes in net1 reset to group -1 */
/*   ResetNetGroup(net2); // All nodes in net2 reset to group -1 */
/*   part = CreateHeaderGroup(); */
/*   lastg = part; */

/*   p = net1;     // nodes in net1 */
/*   for(i=0; i<S1; i++) { */
/*     p = p->next; */
/*     nlist[i] = p; */
/*     glist[i] = CreateGroup(lastg, i); */
/*     lastg = glist[i]; */
/*     AddNodeToGroup(glist[i], p); */
/*   } */
/*   p = net2;     // nodes in net2 */
/*   for(i=0; i<S2; i++) { */
/*     p = p->next; */
/*     nlist[i+S1] = p; */
/*     glist[i+S1] = CreateGroup(lastg, i+S1); */
/*     lastg = glist[i+S1]; */
/*     AddNodeToGroup(glist[i+S1], p); */
/*   } */

/*   // Number of iterations at each temperature */
/*   cicle1 = floor(fac * (double)(nnod *nnod)); */
/*   if (cicle1 < 10) */
/*     cicle1 = 10; */
/*   cicle2 = floor(fac * (double)nnod); */
/*   if (cicle2 < 2) */
/*     cicle2 = 2; */

/*   // START THE SIMULATED ANNEALING */
/*   T = Ti; */
/*   energy = BipartCoModularity(binet, part); */

/*   while ((T > Tf) && (count < limit)) { */

/* /\*     printf("%g %lf %g\n", 1.0/T, energy, T); *\/ */
/*     printf("%g %lf %lf %g %d\n", 1.0/T, energy, */
/* 	   BipartCoModularity(binet, part), T, */
/* 	   CountNonEmptyGroups(part)); */
    
/*     /\* */
/*       Do cicle2 collective change iterations */
/*     *\/ */
/*     if (merge == 1) { */
/*       for (i=0; i<cicle2; i++) { */
	
/* 	// Merge ------------------------------ */
/* 	target = floor(gsl_rng_uniform(gen) * nnod); */
/* 	g1 = nlist[target]->inGroup; */

/* 	if (glist[g1]->size < nnod) { */
/* 	  do { */
/* 	    target = floor(gsl_rng_uniform(gen) * nnod); */
/* 	    g2 = nlist[target]->inGroup; */
/* 	  } while (g1 == g2); */
	  
/* 	  // Calculate dE */
/* 	  dE = 0.0; */
/* 	  nod = glist[g1]->nodeList; */
/* 	  while ((nod = nod->next) != NULL) { */
/* 	    nod2 = glist[g2]->nodeList; */
/* 	    while ((nod2 = nod2->next) != NULL) { */
/* 	      if (nod->ref->ivar1 != nod2->ref->ivar1) { */
/* 		t1 = NodeDegree(nod->ref); */
/* 		t2 = NodeDegree(nod2->ref); */
/* 		dE += (IsThereLink(nod->ref, nod2->node) -  */
/* 		       (double)(t1 * t2) / (double)(L)) /  */
/* 		  (double)L; */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	  // Accept/reject change */
/* 	  if ((dE > 0) || (gsl_rng_uniform(gen) < exp(dE/T))) { */
/* /\* 	    printf("\taccepting merge\n"); *\/ */
/* 	    MergeGroups(glist[g1], glist[g2]); */
/* 	    energy += dE; */
/* 	  } */
/* /\* 	  else { *\/ */
/* /\* 	    printf("\trejecting merge\n"); *\/ */
/* /\* 	  } *\/ */
/* 	} // End of merge move */
	
/* 	// Split ------------------------------ */
/* 	// Look for an empty group */
/* 	g = part; */
/* 	empty = -1; */
/* 	while (((g = g->next) != NULL) && (empty < 0)) { */
/* 	  if (g->size == 0) { */
/* 	    empty = g->label; */
/* 	    break; */
/* 	  } */
/* 	} */

/* 	if (empty >= 0 ) { // if there are no empty groups, do nothing */
/* 	  // Select group */
/* 	  do { */
/* 	    target = floor(gsl_rng_uniform(gen) * (double)nnod); // targ */
/* 							       // node */
/* 	    target = nlist[target]->inGroup;    // target group */
/* 	  } while (glist[target]->size == 1); */

/* 	  // Split the group */
/* 	  if (prob < 0.) */
/* 	    ThermalBipartworkCoSplit(glist[target], glist[empty], */
/* 				    Ti, T, gen); */
/* 	  else */
/* 	    ThermalPercBipartworkCoSplit(glist[target], glist[empty], */
/* 					prob, Ti, T, */
/* 					gen); */
	  
/* 	  // Calculate dE for remerging the groups */
/* 	  dE = 0.0; */
/* 	  nod = glist[target]->nodeList; */
/* 	  while ((nod = nod->next) != NULL) { */
/* 	    nod2 = glist[empty]->nodeList; */
/* 	    while ((nod2 = nod2->next) != NULL) { */
/* 	      if (nod->ref->ivar1 != nod2->ref->ivar1) { */
/* 		t1 = NodeDegree(nod->ref); */
/* 		t2 = NodeDegree(nod2->ref); */
/* 		dE += (IsThereLink(nod->ref, nod2->node) - */
/* 		       (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	      } */
/* 	    } */
/* 	  } */

/* 	  // Accept the change according to "inverse" Metroppolis. */
/* 	  // Inverse means that the algor is applied to the split and */
/* 	  // NOT to the merge! */
/* 	  if ((dE > 0.0) && (gsl_rng_uniform(gen) > exp(-dE/T))) { */
/* 	    // Undo the split */
/* 	    MergeGroups(glist[target],glist[empty]); */
/* /\* 	    printf("\trejecting split\n"); *\/ */
/* 	  } */
/* 	  else{ */
/* 	    // Update energy */
/* /\* 	    printf("\taccepting split\n"); *\/ */
/* 	    energy -= dE; */
/* 	  } */
/* 	} // End of split move */
	
/*       } // End of cicle2 loop */
/*     } // End of 'if merge==1' loop */

/*     /\* */
/*       Do cicle1 individual change iterations */
/*     *\/ */
/*     for (i=0; i<cicle1; i++) { */

/*       // Propose an individual change // */
/*       target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do { */
/* 	newg = floor(gsl_rng_uniform(gen) * nnod); */
/*       } while (newg == oldg); */

/*       // Calculate the change of energy */
/*       dE = 0.0; */
/*       t1 = NodeDegree(nlist[target]); */
/*       // Old group contribution */
/*       nod = glist[oldg]->nodeList; */
/*       while ((nod = nod->next) != NULL) { */
/* 	if (nlist[target]->ivar1 != nod->ref->ivar1) { */
/* 	  t2 = NodeDegree(nod->ref); */
/* 	  dE -= (IsThereLink(nlist[target], nod->node) -  */
/* 		 (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	} */
/*       } */

/*       // New group contribution */
/*       nod = glist[newg]->nodeList; */
/*       while ((nod = nod->next) != NULL) { */
/* 	if (nlist[target]->ivar1 != nod->ref->ivar1) { */
/* 	  t2 = NodeDegree(nod->ref); */
/* 	  dE += (IsThereLink(nlist[target], nod->node) -  */
/* 		 (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	} */
/*       } */

/*       // Accept/reject movement according to Metropolis */
/*       if ((dE > 0) || (gsl_rng_uniform(gen) < exp(dE/T))) { */
/* 	energy += dE; */
/* 	MoveNode(nlist[target],glist[oldg],glist[newg]); */
/*       } */
/*     } */

/*     // Update the no-change counter */
/*     if (energy == energyant) { */
/*       count++; */
      
/*       // If the program is ready to stop (count==limit) but the */
/*       // current partition is not the best one so far, replace the */
/*       // current partition by the best one and continue from there. */
/*       if ((count == limit) && (energy < best_E)) { */
/* 	printf("# Resetting partition\n"); */
/* 	RemovePartition(part); */
/* 	g = part = CopyPartition(best_part); */
/* 	while ((g = g->next) != NULL) { */
/* 	  glist[g->label] = g; */
/* 	  nod = g->nodeList; */
/* 	  while ((nod = nod->next) != NULL) { */
/* 	    nod->ref->inGroup = g->label; */
/* 	  } */
/* 	} */

/* 	// NEED TO MAP binet TO part???? */

/* 	energy = best_E; */
/* 	count = 0; */
/*       } */
/*     } */
/*     else { */
/*       energyant = energy; */
/*       count = 0; */
/*     } */

/*     // Compare the current partition to the best partition so far and */
/*     // save the current if it is better than the best so far. */
/*     if ( energy > best_E ) { */
/*       if ( best_part != NULL ) */
/* 	RemovePartition(best_part); */
/*       // Save the best partition */
/*       best_part = CopyPartition(part); */
/*       // NEED TO MAP binet TO part???? */
/*       best_E = energy; */
/*     } */

/*     // Update the temperature */
/*     T = T * Ts; */

/*   } // End of temperature loop */

/*   RemovePartition(best_part); */
/*   return CompressPart(part); */
/* } */


