/*
  modules.c
  $LastChangedDate$
  $Revision$
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include "prng.h"

#include "tools.h"
#include "graph.h"
#include "modules.h"

#include "bipartite.h"

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
CreateBinet()
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

  /* Initialize the subnetworks */
  last1 = root1 = CreateHeaderGraph();
  last2 = root2 = CreateHeaderGraph();

  /* Go through the file */
  while (!feof(inFile)) {

    /* Read the labels (and weight, if necessary) */
    if (weight_sw == 0) {
      fscanf(inFile,"%s %s\n", label1, label2);
      weight = 1.;
    }
    else {
      fscanf(inFile,"%s %s %lf\n", label1, label2, &weight);
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
      FreeNodeTree(nodeTree, preorder, 0);
    n1 = ntree1->ref;

    nodeTree = CreateNodeTree();
    strcpy(nodeTree->label, label2);
    ntree2 = *(struct node_tree **)tsearch((void *)nodeTree,
					   &dict2,
					   NodeTreeLabelCompare);
    if (ntree2->ref == NULL)
      ntree2->ref = last2 = CreateNodeGraph(last2, label2);
    else
      FreeNodeTree(nodeTree, preorder, 0);
    n2 = ntree2->ref;

    /* Create the links */
    AddAdjacency(n1, n2, 0, add_weight_sw, weight, 0);
    AddAdjacency(n2, n1, 0, add_weight_sw, weight, 0);
  }

  /* Free memory */
  tdestroy(dict1, FreeNodeTree);
  tdestroy(dict2, FreeNodeTree);

  /* Create the bipartite network and return */
  net = CreateBinet();
  net->net1 = root1;
  net->net2 = root2;
  return net;
}

/*
  ---------------------------------------------------------------------
  Free the memory allocated to a bipartite network
  ---------------------------------------------------------------------
*/
void
RemoveBinet(struct binet *net)
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
  ---------------------------------------------------------------------
  Network operations
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
/*
  ---------------------------------------------------------------------
  Invert the networks in binet
  ---------------------------------------------------------------------
*/
struct binet *
InvertBinet(struct binet *net)
{
  struct node_gra *temp;

  temp = net->net1;
  net->net1 = net->net2;
  net->net2 = temp;

  return net;
}

/*
  ---------------------------------------------------------------------
  Create a copy of a bipartite network
  ---------------------------------------------------------------------
*/
struct binet *
CopyBinet(struct binet *binet)
{
  struct binet *copy = NULL;
  struct node_gra *p1 = NULL, *p2 = NULL;
  struct node_gra *last1 = NULL, *last2 = NULL;
  struct node_lis *l = NULL;
  void *dict1 = NULL, *dict2 = NULL;
  struct node_tree *ntree1=NULL, *ntree2=NULL;

  /* Create the copy binet */
  copy = CreateBinet();
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
  tdestroy(dict1, FreeNodeTree);
  tdestroy(dict2, FreeNodeTree);
  
  /* Done */
  return copy;
}

/*
  ---------------------------------------------------------------------
  Project a binet into a weighted one mode network. The network is
  projected into the space of net1, and the weights in the projection
  are the number of common links in the binet.
  ---------------------------------------------------------------------
*/
struct node_gra *
ProjectBinet(struct binet *binet)
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
  tdestroy(dict, FreeNodeTree);
  return projnet;
}




/* /\* */
/*   --------------------------------------------------------------------- */
/*   Build a modular bipartine network. */

/*   mod_size: vector with the sizes of the modules (if NULL, then the */
/*   modules are all of the same size nodpermod). */

/*   nodpermod: number of nodes per module (NOTE: only used if */
/*   mod_size==NULL) */

/*   nmod: number of modules */

/*   S2: number of 'teams' */

/*   mmin (mmax): minimum (maximum) number of 'agents' per 'team'. Team */
/*   size is determined by drawing a randomly generated integer uniformly */
/*   distributed between mmin and mmax. */

/*   p: probability of selecting agents of a given color for a team */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* struct binet *BuildModularBipartiteNetwork(int *mod_size, */
/* 					   int nodpermod, */
/* 					   int nmod, */
/* 					   double *col_prob, */
/* 					   int S2, */
/* 					   int mmin, int mmax, */
/* 					   double geom_p, */
/* 					   double p, */
/* 					   struct prng *gen) */
/* { */
/*   int i, j; */
/*   struct binet *net = NULL; */
/*   struct node_gra *root1 = NULL; */
/*   struct node_gra *root2 = NULL; */
/*   struct node_gra *list1[maxim_int]; */
/*   struct node_gra *list2[maxim_int]; */
/*   struct node_gra *last_add1 = NULL; */
/*   struct node_gra *last_add2 = NULL; */
/*   int orig, dest; */
/*   int teamcol, target; */
/*   int S1 = 0; */
/*   int *module, *mod_starting, n_team_col; */
/*   int m; */
/*   int free_mod_size = 0; */
/*   double dice, cum; */

/*   // Define mod_size if not defined yet */
/*   if (mod_size == NULL) { */
/*     free_mod_size = 1; // Remember to free memory at the end */
/*     mod_size = allocate_i_vec(nmod); */
/*     for (i=0; i<nmod; i++) { */
/*       mod_size[i] = nodpermod; */
/*     } */
/*   } */

/*   // Calculate the number of nodes S1 and allocate memory */
/*   S1 = 0; */
/*   for (i=0; i<nmod; i++) { */
/*     S1 += mod_size[i]; */
/*   } */
/*   module = allocate_i_vec(S1); */

/*   // Assign each node to a module */
/*   S1 = 0; */
/*   mod_starting = allocate_i_vec(nmod); */
/*   for (i=0; i<nmod; i++) { */
/*     mod_starting[i] = S1;  // Number of the first node in the module */
/*     for (j=0; j<mod_size[i]; j++) { */
/*       module[S1++] = i; */
/*     } */
/*   } */

/*   // Create the two networks */
/*   last_add1 = root1 = CreateHeaderGraph(); */
/*   last_add2 = root2 = CreateHeaderGraph(); */

/*   for (i = 0; i<S1; i++) { */
/*     list1[i] = CreateNodeGraph(last_add1, i, i, i); */
/*     list1[i]->inGroup = module[i]; */
/*     list1[i]->ivar1 = 0; */
/*     last_add1 = list1[i]; */
/*   } */
/*   for (i = 0; i<S2; i++) { */
/*     list2[i] = CreateNodeGraph(last_add2, i, i, i); */
/*     list2[i]->ivar1 = 1; */
/*     last_add2 = list2[i]; */
/*   } */

/*   // Create the links */
/*   for (i = 0; i<S2; i++) {  // loop over teams */

/*     // Determine randomly which is the "color" of the team (according */
/*     // to the probabilities in col_prob if available, or with equal */
/*     // probability for all colors otherwise), that is, the module that */
/*     // will contribute, in principle, more nodes. */
/*     if (col_prob != NULL) { */
/*       dice = prng_get_next(gen); */
/*       cum = 0.0; */
/*       teamcol = -1; */
/*       while (cum < dice) { */
/* 	teamcol++; */
/* 	cum += col_prob[teamcol]; */
/*       } */
/*     } */
/*     else { */
/*       teamcol = floor(prng_get_next(gen) * nmod); */
/*     } */
/*     list2[i]->inGroup = teamcol; */
/*     n_team_col = 0; // No nodes of the team color added so far */

/*     // Determine the number of actors in the team */
/*     if (mmin < 0) { // use a geometric distribution */
/*       m = geometric_dist_val(geom_p, gen); */
/*     } */
/*     else if (mmax != mmin) */
/*       m = floor(prng_get_next(gen) * (double)(mmax+1-mmin) + mmin); */
/*     else */
/*       m = mmin; */

/*     for (j=0; j<m; j++) {  // loop over spots in a team */

/*       if (prng_get_next(gen) < p &&  */
/* 	  n_team_col < mod_size[teamcol]) { // select node with the */
/* 					    // module's color */
/* 	n_team_col++; */
/* 	do { */
/* 	  target = mod_starting[teamcol] + floor(prng_get_next(gen) * */
/* 						 mod_size[teamcol]); */
/* 	} while (ExistLink(list1[target], list2[i]) == 1); */
/*       } */
/*       else {                 // choose node at random */
/* 	do { */
/* 	  target = floor(prng_get_next(gen) * S1); */
/* 	} while (ExistLink(list1[target], list2[i]) == 1); */
/* 	if (module[target] == teamcol) { */
/* 	  n_team_col++; */
/* 	} */
/*       } */
/*       // Make the link */
/*       AddAdjacencyFull(list1[target],list2[i],10); */
/*       AddAdjacencyFull(list2[i],list1[target],10); */

/*     } // end of loop over spots */
      
/*   } // end of loop over teams */
  
/*   // Create the bipartite network */
/*   net = CreateBinet(); */
/*   net->net1 = root1; */
/*   net->net2 = root2; */

/*   // Free memory */
/*   free_i_vec(module); */
/*   free_i_vec(mod_starting); */
/*   if (free_mod_size == 1) { */
/*     free_i_vec(mod_size); */
/*     mod_size = NULL; */
/*   } */

/*   // Done */
/*   return net; */
/* } */


/* /\* */
/*   --------------------------------------------------------------------- */
/*   Counts the links in a bipartite network */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* int TotalNLinksBinet(struct binet *binet) */
/* { */
/*   struct node_gra *p = binet->net1; */
/*   int nlink = 0; */

/*   while ((p = p->next) !=  NULL) */
/*     nlink +=  CountLinks(p); */

/*   return nlink; */
/* } */


/* /\* */
/*   --------------------------------------------------------------------- */
/*   Randomize links in a bipartite network */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* struct binet *RandomizeBinet(struct binet *binet, */
/* 			     double times, struct prng *gen) */
/* { */
/*   int i; */
/*   int nlink, niter, coun = 0; */
/*   int target1, target2; */
/*   struct node_gra *p = NULL; */
/*   struct node_lis *l = NULL; */
/*   struct node_gra *n1, *n2, *n3, *n4; */
/*   struct node_gra **ori, **des; */

/*   // Build the link lists (one for link origins and one for ends) */
/*   nlink = TotalNLinksBinet(binet); */
/*   niter = ceil(times * (double)nlink); */

/*   ori = (struct node_gra **)calloc(nlink,sizeof(struct node_gra *)); */
/*   des = (struct node_gra **)calloc(nlink,sizeof(struct node_gra *)); */

/*   p = binet->net1; */
/*   while ((p = p->next) !=  NULL) { */
/*     l = p->neig; */
/*     while ((l = l->next) !=  NULL) { */
/*       ori[coun] = p; */
/*       des[coun] = l->ref; */
/*       coun++; */
/*     } */
/*   } */

/*   if(coun !=  nlink) */
/*     printf("Error in RandomizeBinet: coun !=  nlink!!\n"); */

/*   // Randomize the network */
/*   for (i=0; i<niter; i++) { */
    
/*     // select the two links (four nodes) to swap */
/*     do { */
/*       target1 = floor(prng_get_next(gen) * (double)nlink); */
/*       n1 = ori[target1]; */
/*       n2 = des[target1]; */

/*       do { */
/* 	target2 = floor(prng_get_next(gen) * (double)nlink); */
/* 	n3 = ori[target2]; */
/* 	n4 = des[target2]; */
/*       } while (n1 == n3 || n2 == n4); */

/*     } while (IsThereLink(n1, n4->num) == 1 || */
/* 	     IsThereLink(n2, n3->num) == 1); */

/*     // switch the link */
/*     RemoveLink(n1, n2); */
/*     RemoveLink(n3, n4); */
/*     AddAdjacencyFull(n1, n4, 1); */
/*     AddAdjacencyFull(n4, n1, 1); */
/*     AddAdjacencyFull(n3, n2, 1); */
/*     AddAdjacencyFull(n2, n3, 1); */

/*     ori[target1] = n1; */
/*     des[target1] = n4; */
/*     ori[target2] = n3; */
/*     des[target2] = n2; */
/*   } */

/*   free(ori); */
/*   free(des); */

/*   return binet; */
/* } */


/* /\* */
/*   --------------------------------------------------------------------- */
/*   Print a bipartite network in Pajek format */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* PrintPajekBinet(struct binet *net) */
/* { */
/*   struct node_gra *p = NULL; */
/*   struct node_lis *li; */
/*   FILE *fit; */
/*   int trans1[maxim_int]; */
/*   int trans2[maxim_int]; */
/*   int i,co; */
/*   int Stot, S1, S2; */

/*   S1 = CountNodes(net->net1); */
/*   S2 = CountNodes(net->net2); */
/*   Stot  = S1 + S2; */

/*   for(i = 0;i<maxim_int;i++){ */
/*     trans1[i] = 0; */
/*     trans2[i] = 0; */
/*   } */

/*   fit = fopen("bigraph.net","w"); */

/*   fprintf(fit,"*Vertices %d %d\n", Stot, S1); */

/*   co = 1; */

/*   // Print the nodes */
/*   p = net->net1; */
/*   while(p->next != NULL){ */
/*     p = p->next; */

/*     trans1[p->num] = co++; */
/*     fprintf(fit,"%d \"A%d\"\n",trans1[p->num],p->num+1); */
/*   } */

/*   p = net->net2; */
/*   while(p->next != NULL){ */
/*     p = p->next; */

/*     trans2[p->num] = co++; */
/*     fprintf(fit,"%d \"B%d\"\n",trans2[p->num],p->num+1); */
/*   } */

/*   // Print the links */
/*   fprintf(fit,"*Edges\n"); */
/*   p = net->net1; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     li = p->neig; */
/*     while(li->next != NULL){ */
/*       li = li->next; */
/*       fprintf(fit,"%d %d\n",trans1[p->num],trans2[li->node]); */
/*     } */
/*   } */

/*   fclose(fit); */
/* } */


/* /\* */
/*   --------------------------------------------------------------------- */
/*   Print a bipartite network in Pajek format */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* PrintPajekBinetFName(struct binet *net, char *fname) */
/* { */
/*   struct node_gra *p = NULL; */
/*   struct node_lis *li; */
/*   FILE *fit; */
/*   int trans1[maxim_int]; */
/*   int trans2[maxim_int]; */
/*   int i,co; */
/*   int Stot, S1, S2; */

/*   S1 = CountNodes(net->net1); */
/*   S2 = CountNodes(net->net2); */
/*   Stot  = S1 + S2; */

/*   for(i = 0;i<maxim_int;i++){ */
/*     trans1[i] = 0; */
/*     trans2[i] = 0; */
/*   } */

/*   fit = fopen(fname, "w"); */

/*   fprintf(fit,"*Vertices %d %d\n", Stot, S1); */

/*   co = 1; */

/*   // Print the nodes */
/*   p = net->net1; */
/*   while((p = p->next) != NULL){ */
/*     trans1[p->num] = co++; */
/*     fprintf(fit,"%d \"A%d\"           %lf %lf %lf\n", */
/* 	    trans1[p->num], p->num+1, */
/* 	    p->coorX, p->coorY, 0.5); */
/*   } */
/*   p = net->net2; */
/*   while((p = p->next) != NULL){ */
/*     trans2[p->num] = co++; */
/*     fprintf(fit,"%d \"B%d\"           %lf %lf %lf\n", */
/* 	    trans2[p->num], p->num+1, */
/* 	    p->coorX, p->coorY, 0.5); */
/*   } */

/*   // Print the links */
/*   fprintf(fit,"*Arcs\n"); */
/*   fprintf(fit,"*Edges\n"); */
/*   p = net->net1; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     li = p->neig; */
/*     while(li->next != NULL){ */
/*       li = li->next; */
/*       fprintf(fit,"%d %d\n",trans1[p->num],trans2[li->node]); */
/*     } */
/*   } */

/*   fclose(fit); */
/* } */




/* /\* */
/*   --------------------------------------------------------------------- */
/*   Calculate the modularity of a bipartite network */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* double BinetModularity(struct binet *binet, struct group *part) */
/* { */
/*   struct node_gra *p = binet->net2; */
/*   struct node_lis *n1, *n2; */
/*   double bimod = 0.0; */
/*   int c12; */
/*   int t1, t2; */
/*   double sms, sms2; */

/*   sms = sms2 = 0.0; */

/*   // Calculate sms=sum(m_s) and sms2=sum(m_s^2) */
/*   while (p->next != NULL) { */
/*     p = p->next; */
/*     sms += (double)CountLinks(p); */
/*     sms2 += (double)(CountLinks(p) * CountLinks(p)); */
/*   } */

/*   // Calculate the modularity */
/*   while (part->next != NULL){ */
/*     part = part->next; */

/*     n1 = part->nodeList; */
/*     while ((n2 = n1 = n1->next) != NULL) { */
/*       while ((n2 = n2->next) != NULL) { */
/*   	c12 = NCommonLinksBipart(n1->ref, n2->ref); */
/* 	t1 = CountLinks(n1->ref); */
/* 	t2 = CountLinks(n2->ref); */

/* 	bimod += (double)c12 / (sms2 - sms) - */
/* 	  (double)(t1 * t2) / (sms * sms); */
/*       } */
/*     } */
/*   } */

/*   bimod *= 2.; */

/*   return bimod; */
/* } */


/* /\* */
/*   --------------------------------------------------------------------- */
/*   Calculate the modularity of a co-partition of both subnetworks of a */
/*   bipartite network */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* double BinetCoModularity(struct binet *binet, struct group *part) */
/* { */
/*   struct node_lis *n1, *n2; */
/*   double bimod = 0.0; */
/*   int t1, t2; */
/*   int L; */
/* /\*   int S1, S2; *\/ */

/* /\*   S1 = CountNodes(binet->net1); *\/ */
/* /\*   S2 = CountNodes(binet->net2); *\/ */
/*   L = TotalNLinksBinet(binet); */

/*   // Calculate the modularity */
/*   while ((part = part->next) != NULL){ */
/*     n1 = part->nodeList; */
/*     while ((n2 = n1 = n1->next) != NULL) { */
/*       while ((n2 = n2->next) != NULL) { */
/* 	if ( n1->ref->ivar1 != n2->ref->ivar1 ) { */
/* 	  t1 = CountLinks(n1->ref); */
/* 	  t2 = CountLinks(n2->ref); */
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
/*   Splits one group into two new groups (thermalized with the outer */
/*   SA), checking first if the network contains disconnected clusters. */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* ThermalPercBinetworkSplit(struct group *target_g, */
/* 			  struct group *empty_g, */
/* 			  double prob, double Ti, double Tf, */
/* 			  double **cmat, double msfac, */
/* 			  struct prng *gen) */
/* { */
/*   struct group *glist[2], *g = NULL, *split = NULL; */
/*   struct node_gra *nlist[maxim_int]; */
/*   int trans[maxim_int]; */
/*   struct node_lis *p = NULL, *lastp = NULL; */
/*   int nnod = 0; */
/*   int i; */
/*   int n1, n2, t1, t2; */
/*   int target, oldg, newg; */
/*   double dE = 0.0; */
/*   double T, Ts = 0.95; */
/*   double dice; */
/*   struct node_gra *net; */
/*   struct binet *binet; */
  
/*   /\* */
/*     Map the groups and the nodes */
/*   *\/ */
/*   glist[0] = target_g; */
/*   glist[1] = empty_g; */

/*   p = target_g->nodeList; */
/*   while ((p = p->next) != NULL) { */
/*     nlist[nnod] = p->ref; */
/*     trans[p->node] = nnod++; */
/*   } */

/*   /\* */
/*     Check if the network is disconnected */
/*   *\/ */
/*   // Build a network from the nodes in the target group */
/*   binet = CreateBinet(); */
/*   binet->net1 = BuildNetFromPart(target_g); */
/*   binet->net2 = CreateHeaderGraph(); */
/*   net = ProjectBinet(binet); */
/*   split = ClustersPartition(net); */

/*   /\*  */
/*      If the network is disconnected */
/*   *\/ */
/* /\*   printf("Clusters: %d\n", CountGroups(split)); *\/ */
/*   if (CountGroups(split) > 1 && prng_get_next(gen) < prob) { */
/* /\*     printf("Going perc\n"); *\/ */
/*     // Move the nodes in some (half) clusters to empty module */
/*     g = split; */
/*     while ((g = g->next) != NULL) { */
/*       if (prng_get_next(gen) < 0.5) { // With prob=0.5 move to empty */
/* 	p = g->nodeList; */
/* 	while ((p = p->next) != NULL) { */
/* 	  MoveNode(nlist[trans[p->node]], target_g, empty_g); */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   /\* */
/*     Else if the network is connected */
/*   *\/ */
/*   else { */
/* /\*     printf("Going SA\n"); *\/ */
/*     // Randomly assign the nodes to the groups */
/*     lastp = p = target_g->nodeList; */
/*     while ((p = p->next) != NULL) { */
/*       dice = prng_get_next(gen); */
/*       if (dice < 0.5) { */
/* 	MoveNode(p->ref, target_g, empty_g); */
/* 	p = lastp; */
/*       } */
/*       else { */
/* 	lastp = p; */
/*       } */
/*     } */

/*     // Do SA to "optimize" the splitting */
/*     T = Ti; */
/*     while (T >= Tf) { */

/*       for (i=0; i<nnod; i++) { */
/* 	target = floor(prng_get_next(gen) * (double)nnod); */
/* 	if (nlist[target]->inGroup == target_g->label) */
/* 	  oldg = 0; */
/* 	else */
/* 	  oldg = 1; */
/* 	newg = 1 - oldg; */
	
/* 	// Calculate the change of energy */
/* 	dE = 0.0; */
/* 	n1 = nlist[target]->num; */
/* 	t1 = CountLinks(nlist[target]); */
	
/* 	// Old group */
/* 	p = glist[oldg]->nodeList; */
/* 	while ((p = p->next) != NULL) { */
/* 	  n2 = p->node; */
/* 	  if (n2 != n1) { */
/* 	    t2 = CountLinks(p->ref); */
/* 	    dE -= 2. * (cmat[n1][n2] - t1 * t2 * msfac); */
/* 	  } */
/* 	} */

/* 	// New group */
/* 	p = glist[newg]->nodeList; */
/* 	while ((p = p->next) != NULL) { */
/* 	  n2 = p->node; */
/* 	  t2 = CountLinks(p->ref); */
/* 	  dE += 2. * (cmat[n1][n2] - t1 * t2 * msfac); */
/* 	} */
	
/* 	// Accept the change according to the Boltzman factor */
/* 	if( (dE >= 0.0) || (prng_get_next(gen) < exp(dE/T)) ){ */
/* 	  MoveNode(nlist[target], glist[oldg], glist[newg]); */
/* 	} */
/*       } */
      
/*       T = T * Ts; */
/*     } // End of temperature loop */

/*   } // End of else (network is connected) */
  
/*   RemovePartition(split); */
/*   RemoveGraph(net); */
/*   RemoveBinet(binet); */
/* } */


/* /\* */
/*   --------------------------------------------------------------------- */
/*   Splits one group into two new groups (thermalized with the outer SA) */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* ThermalBinetworkSplit(struct group *target_g, struct group *empty_g, */
/* 		      double Ti, double Tf, */
/* 		      double **cmat, double msfac, */
/* 		      struct prng *gen) */
/* { */
/*   struct group *glist[2]; */
/*   struct node_gra *nlist[maxim_int]; */
/*   struct node_lis *p = NULL, *lastp = NULL; */
/*   int nnod = 0; */
/*   int i; */
/*   int n1, n2, t1, t2; */
/*   int target, oldg, newg; */
/*   double dE = 0.0; */
/*   double T, Ts = 0.95; */
/*   double dice; */

/*   // Map the groups */
/*   glist[0] = target_g; */
/*   glist[1] = empty_g; */

/*   // Randomly assign the nodes to the groups */
/*   lastp = p = target_g->nodeList; */
/*   while ((p = p->next) != NULL) { */
/*     nlist[nnod++] = p->ref; */
/*     dice = prng_get_next(gen); */
/*     if (dice < 0.5) { */
/*       MoveNode(p->ref, target_g, empty_g); */
/*       p = lastp; */
/*     } */
/*     else { */
/*       lastp = p; */
/*     } */
/*   } */

/*   // Do SA to "optimize" the splitting */
/*   T = Ti; */
/*   while (T >= Tf) { */

/*     for (i=0; i<nnod; i++) { */
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       if (nlist[target]->inGroup == target_g->label) */
/* 	oldg = 0; */
/*       else */
/* 	oldg = 1; */
/*       newg = 1 - oldg; */
	
/*       // Calculate the change of energy */
/*       dE = 0.0; */
/*       n1 = nlist[target]->num; */
/*       t1 = CountLinks(nlist[target]); */
	
/*       // Old group */
/*       p = glist[oldg]->nodeList; */
/*       while ((p = p->next) != NULL) { */
/* 	n2 = p->node; */
/* 	if (n2 != n1) { */
/* 	  t2 = CountLinks(p->ref); */
/* 	  dE -= 2. * (cmat[n1][n2] - t1 * t2 * msfac); */
/* 	} */
/*       } */

/*       // New group */
/*       p = glist[newg]->nodeList; */
/*       while ((p = p->next) != NULL) { */
/* 	n2 = p->node; */
/* 	t2 = CountLinks(p->ref); */
/* 	dE += 2. * (cmat[n1][n2] - t1 * t2 * msfac); */
/*       } */

/*       // Accept the change according to the Boltzman factor */
/*       if( (dE >= 0.0) || (prng_get_next(gen) < exp(dE/T)) ){ */
/* 	MoveNode(nlist[target], glist[oldg], glist[newg]); */
/*       } */
/*     } */

/*     T = T * Ts; */
/*   } // End of temperature loop */
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
/*   binet = CreateBinet(); */
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
/* ThermalBinetworkCoSplit(struct group *target_g, struct group *empty_g, */
/* 			double Ti, double Tf, */
/* 			struct prng *gen) */
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
/*   L = TotalNLinksBinet(module_binet); */
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
/*     nlinks[p->trans] = CountLinks(p); */
/*     dice = prng_get_next(gen); */
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
/*     nlinks[p->trans] = CountLinks(p); */
/*     dice = prng_get_next(gen); */
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
/*       target = floor(prng_get_next(gen) * (double)nnod); */
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
/*       if( (dE >= 0.0) || (prng_get_next(gen) < exp(dE/T)) ){ */
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
/*   RemoveBinet(module_binet); */
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
/* struct group *BinetClustersPartition(struct binet *binet) */
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
/*   Same as ThermalBinetCoSplit, but checks first if there are */
/*   disconnected components. If indeed the module is disconnected, */
/*   proposes a split using this fact. */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* ThermalPercBinetworkCoSplit(struct group *target_g, struct group *empty_g, */
/* 			    double prob, double Ti, double Tf, */
/* 			    struct prng *gen) */
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
/*   L = TotalNLinksBinet(module_binet); */
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
/*   temp_part = BinetClustersPartition(module_binet); */

/*   if (CountGroups(temp_part) > 1 && prng_get_next(gen) < prob) { */
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
/*       node1 = floor(prng_get_next(gen) * (double)nnod); */
/*       mod1 = nlist[node1]->inGroup; */
/*       do { */
/* 	node2 = floor(prng_get_next(gen) * (double)nnod); */
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
/*       nlinks[p->trans] = CountLinks(p); */
/*       dice = prng_get_next(gen); */
/*       if (dice < 0.5) { */
/* 	AddNodeToGroup(glist[0], p); */
/*       } */
/*       else { */
/* 	AddNodeToGroup(glist[1], p); */
/*       } */
/*     } */
/*     p = module_binet->net2; */
/*     while ((p = p->next) != NULL) { */
/*       nlinks[p->trans] = CountLinks(p); */
/*       dice = prng_get_next(gen); */
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
/* 	target = floor(prng_get_next(gen) * (double)nnod); */
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
/* 	if( (dE >= 0.0) || (prng_get_next(gen) < exp(dE/T)) ){ */
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
/*   RemoveBinet(module_binet); */
/*   RemovePartition(temp_part); */
/*   free_i_vec(nlinks); */
/* } */


/* /\* */
/*   --------------------------------------------------------------------- */
/*   Identify modules in the net1 network of a bipartite network using */
/*   simulated annealing. */

/*   Ti: Initial temperature for the SA */
/*   Tf: Final temperature for the SA */
/*   Ts: Cooling factor */
/*   fac: Iteration factor */
/*   merge: implement collective moves (merge=1) or not (merge=0) */
/*   prob: only used if merge==1. Use percolation split with probability */
/*   prob (prob>0) or do not use percolation (prob<0) */
/*   --------------------------------------------------------------------- */
/* *\/ */
/* struct group *SACommunityIdentBinet(struct binet *binet, */
/* 				    double Ti, double Tf, double Ts, */
/* 				    double fac, int merge, double prob, */
/* 				    struct prng *gen) */
/* { */
/*   int i; */
/*   struct node_gra *net1 = binet->net1; */
/*   struct node_gra *net2 = binet->net2; */
/*   double sms = 0.0, sms2 = 0.0, msfac; */
/*   struct node_gra *p, *p2; */
/*   int nnod; */
/*   struct group *part = NULL, *g = NULL, *split = NULL; */
/*   struct node_gra *nlist[maxim_int]; */
/*   struct group *glist[maxim_int]; */
/*   int cicle1, cicle2; */
/*   int count = 0, limit = 25; */
/*   double energy, energyant, dE, e; */
/*   double T; */
/*   int target, oldg, newg; */
/*   double **cmat; */
/*   int g1, g2, empty; */
/*   struct node_lis *nod, *nod2; */
/*   double t1, t2; */
/*   struct group *best_part = NULL; */
/*   double best_E = -100.0; */

/*   cmat = allocate_d_mat(maxim_int, maxim_int); */

/*   // Calculate s2ms=(sum m_s)^2 and sms2=sum(m_s^2) */
/*   p = net2; */
/*   while ((p = p->next) != NULL) { */
/*     sms += (double)CountLinks(p); */
/*     sms2 += (double)(CountLinks(p) * CountLinks(p)); */
/*   } */
/*   msfac = 1. / (sms * sms); */

/*   // Create the groups and assign each node to one group */
/*   nnod = CountNodes(net1); */
/*   part = CreateHeaderGroup(); */
/*   ResetNetGroup(net1); // All nodes reset to group -1 */
/*   p = net1->next; */

/*   nlist[0] = p; */
/*   glist[0] = CreateGroup(part, 0); */
/*   AddNodeToGroup(glist[0], p); */
  
/*   for(i=1; i<nnod; i++) { */
/*     p = p->next; */
/*     nlist[i] = p; */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*     AddNodeToGroup(glist[i],p); */
/*   } */

/*   // Calculate the c matrix of concurrences */
/*   p = net1; */
/*   while(p->next != NULL){ */
/*     p = p->next; */

/*     p2 = net1; */
/*     while(p2->next != NULL){ */
/*       p2 = p2->next; */

/*       cmat[p->num][p2->num] = (double)NCommonLinksBipart(p, p2) / */
/* 	(sms2 - sms); */
/*     } */
/*   } */

/*   // Number of iterations at each temperature */
/*   if (fac*(double)(nnod*nnod) < 10) */
/*     cicle1 = 10; */
/*   else */
/*     cicle1 = floor(fac*(double)(nnod*nnod)); */

/*   if (fac*(double)nnod < 2) */
/*     cicle2 = 2; */
/*   else */
/*     cicle2 = floor(fac*(double)nnod); */

/*   // START THE SIMULATED ANNEALING */
/*   T = Ti; */
/*   energy = BinetModularity(binet, part); */

/*   while ((T > Tf) && (count < limit)) { */

/*     printf("%g %lf %g\n", 1.0/T, energy, T); */
/* /\*     printf("%g %lf %lf %g\n", 1.0/T, energy, *\/ */
/* /\* 	   BinetModularity(binet, part), T); *\/ */
    
/*     /\* */
/*       Do cicle2 collective change iterations */
/*     *\/ */
/*     if (merge == 1) { */
/*       for (i=0; i<cicle2; i++){ */
	
/* 	// Merge ------------------------------ */
/* 	target = floor(prng_get_next(gen) * nnod); */
/* 	g1 = nlist[target]->inGroup; */

/* 	if (glist[g1]->size < nnod) { */
/* 	  do { */
/* 	    target = floor(prng_get_next(gen) * nnod); */
/* 	    g2 = nlist[target]->inGroup; */
/* 	  } while (g1 == g2); */
	  
/* 	  // Calculate dE */
/* 	  dE = 0.0; */
/* 	  nod = glist[g1]->nodeList; */
/* 	  while ((nod = nod->next) != NULL) { */
/* 	    nod2 = glist[g2]->nodeList; */
/* 	    while ((nod2 = nod2->next) != NULL) { */
/* 	      t1 = CountLinks(nod->ref); */
/* 	      t2 = CountLinks(nod2->ref); */
/* 	      dE += 2. * (cmat[nod->node][nod2->node] - */
/* 			  t1 * t2 * msfac); */
/* 	    } */
/* 	  } */
/* 	  // Accept/reject change */
/* 	  if ((dE > 0) || (prng_get_next(gen) < exp(dE/T))) { */
/* 	    MergeGroups(glist[g1], glist[g2]); */
/* 	    energy += dE; */
/* 	  } */
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
/* 	    target = floor(prng_get_next(gen) * (double)nnod); // node */
/* 	    target = nlist[target]->inGroup;    // target group */
/* 	  } while (glist[target]->size == 1); */

/* 	  // Split the group */
/* 	  if (prob < 0.) */
/* 	    ThermalBinetworkSplit(glist[target], glist[empty], */
/* 				     Ti, T, cmat, msfac, gen); */
/* 	  else */
/* 	    ThermalPercBinetworkSplit(glist[target], glist[empty], */
/* 				      prob, Ti, T, */
/* 				      cmat, msfac, gen); */
	  
/* 	  // Calculate dE for remerging the groups */
/* 	  dE = 0.0; */
/* 	  nod = glist[target]->nodeList; */
/* 	  while ((nod = nod->next) != NULL) { */
/* 	    nod2 = glist[empty]->nodeList; */
/* 	    while ((nod2 = nod2->next) != NULL) { */
/* 	      t1 = CountLinks(nod->ref); */
/* 	      t2 = CountLinks(nod2->ref); */
/* 	      dE += 2. * (cmat[nod->node][nod2->node] - */
/* 			  t1 * t2 * msfac); */
/* 	    } */
/* 	  } */

/* 	  // Accept the change according to "inverse" Metroppolis. */
/* 	  // Inverse means that the algor is applied to the split and */
/* 	  // NOT to the merge! */
/* 	  if ((dE > 0.0) && (prng_get_next(gen) > exp(-dE/T))) { */
/* 	    // Undo the split */
/* 	    MergeGroups(glist[target],glist[empty]); */
/* 	  } */
/* 	  else{ */
/* 	    // Update energy */
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
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do { */
/* 	newg = floor(prng_get_next(gen) * nnod); */
/*       } while (newg == oldg); */

/*       // Calculate the change of energy */
/*       dE = 0.0; */
/*       t1 = CountLinks(nlist[target]); */
/*       // Old group contribution */
/*       nod = glist[oldg]->nodeList; */
/*       while ((nod = nod->next) != NULL) { */
/* 	t2 = CountLinks(nod->ref); */
/* 	dE -= 2. * (cmat[nlist[target]->num][nod->node] - */
/* 		    t1 * t2 * msfac); */
/*       } */

/*       // New group contribution */
/*       nod = glist[newg]->nodeList; */
/*       while ((nod = nod->next) != NULL) { */
/* 	t2 = CountLinks(nod->ref); */
/* 	dE += 2. * (cmat[nlist[target]->num][nod->node] - */
/* 		    t1 * t2 * msfac); */
/*       } */
/*       dE += 2. * (cmat[nlist[target]->num][nlist[target]->num] - */
/* 		  t1 * t1 * msfac); */

/*       // Accept/reject movement according to Metropolis */
/*       if ((dE > 0) || (prng_get_next(gen) < exp(dE/T))) { */
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
/* 	} */

/* 	MapPartToNet(part, binet->net1); */
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
/*       MapPartToNet(part, binet->net1); // MUST DO this after copying a */
/* 				       // part! */
/*       best_E = energy; */
/*     } */

/*     // Update the temperature */
/*     T = T * Ts; */

/*   } // End of temperature loop */

/*   free_d_mat(cmat, maxim_int); */
/*   RemovePartition(best_part); */
/*   return CompressPart(part); */
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
/* struct group *SACommunityCoIdentBinet(struct binet *binet, */
/* 				      double Ti, double Tf, double Ts, */
/* 				      double fac, */
/* 				      int merge, */
/* 				      double prob, */
/* 				      struct prng *gen) */
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
/*   L = TotalNLinksBinet(binet); */
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
/*   energy = BinetCoModularity(binet, part); */

/*   while ((T > Tf) && (count < limit)) { */

/* /\*     printf("%g %lf %g\n", 1.0/T, energy, T); *\/ */
/*     printf("%g %lf %lf %g %d\n", 1.0/T, energy, */
/* 	   BinetCoModularity(binet, part), T, */
/* 	   CountNonEmptyGroups(part)); */
    
/*     /\* */
/*       Do cicle2 collective change iterations */
/*     *\/ */
/*     if (merge == 1) { */
/*       for (i=0; i<cicle2; i++) { */
	
/* 	// Merge ------------------------------ */
/* 	target = floor(prng_get_next(gen) * nnod); */
/* 	g1 = nlist[target]->inGroup; */

/* 	if (glist[g1]->size < nnod) { */
/* 	  do { */
/* 	    target = floor(prng_get_next(gen) * nnod); */
/* 	    g2 = nlist[target]->inGroup; */
/* 	  } while (g1 == g2); */
	  
/* 	  // Calculate dE */
/* 	  dE = 0.0; */
/* 	  nod = glist[g1]->nodeList; */
/* 	  while ((nod = nod->next) != NULL) { */
/* 	    nod2 = glist[g2]->nodeList; */
/* 	    while ((nod2 = nod2->next) != NULL) { */
/* 	      if (nod->ref->ivar1 != nod2->ref->ivar1) { */
/* 		t1 = CountLinks(nod->ref); */
/* 		t2 = CountLinks(nod2->ref); */
/* 		dE += (IsThereLink(nod->ref, nod2->node) -  */
/* 		       (double)(t1 * t2) / (double)(L)) /  */
/* 		  (double)L; */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	  // Accept/reject change */
/* 	  if ((dE > 0) || (prng_get_next(gen) < exp(dE/T))) { */
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
/* 	    target = floor(prng_get_next(gen) * (double)nnod); // targ */
/* 							       // node */
/* 	    target = nlist[target]->inGroup;    // target group */
/* 	  } while (glist[target]->size == 1); */

/* 	  // Split the group */
/* 	  if (prob < 0.) */
/* 	    ThermalBinetworkCoSplit(glist[target], glist[empty], */
/* 				    Ti, T, gen); */
/* 	  else */
/* 	    ThermalPercBinetworkCoSplit(glist[target], glist[empty], */
/* 					prob, Ti, T, */
/* 					gen); */
	  
/* 	  // Calculate dE for remerging the groups */
/* 	  dE = 0.0; */
/* 	  nod = glist[target]->nodeList; */
/* 	  while ((nod = nod->next) != NULL) { */
/* 	    nod2 = glist[empty]->nodeList; */
/* 	    while ((nod2 = nod2->next) != NULL) { */
/* 	      if (nod->ref->ivar1 != nod2->ref->ivar1) { */
/* 		t1 = CountLinks(nod->ref); */
/* 		t2 = CountLinks(nod2->ref); */
/* 		dE += (IsThereLink(nod->ref, nod2->node) - */
/* 		       (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	      } */
/* 	    } */
/* 	  } */

/* 	  // Accept the change according to "inverse" Metroppolis. */
/* 	  // Inverse means that the algor is applied to the split and */
/* 	  // NOT to the merge! */
/* 	  if ((dE > 0.0) && (prng_get_next(gen) > exp(-dE/T))) { */
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
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do { */
/* 	newg = floor(prng_get_next(gen) * nnod); */
/*       } while (newg == oldg); */

/*       // Calculate the change of energy */
/*       dE = 0.0; */
/*       t1 = CountLinks(nlist[target]); */
/*       // Old group contribution */
/*       nod = glist[oldg]->nodeList; */
/*       while ((nod = nod->next) != NULL) { */
/* 	if (nlist[target]->ivar1 != nod->ref->ivar1) { */
/* 	  t2 = CountLinks(nod->ref); */
/* 	  dE -= (IsThereLink(nlist[target], nod->node) -  */
/* 		 (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	} */
/*       } */

/*       // New group contribution */
/*       nod = glist[newg]->nodeList; */
/*       while ((nod = nod->next) != NULL) { */
/* 	if (nlist[target]->ivar1 != nod->ref->ivar1) { */
/* 	  t2 = CountLinks(nod->ref); */
/* 	  dE += (IsThereLink(nlist[target], nod->node) -  */
/* 		 (double)(t1 * t2) / (double)(L)) / (double)L; */
/* 	} */
/*       } */

/*       // Accept/reject movement according to Metropolis */
/*       if ((dE > 0) || (prng_get_next(gen) < exp(dE/T))) { */
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


