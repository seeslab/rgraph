/*
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

#define EPSILON_MOD 1.e-6

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Group creation and memory allocation
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Create an empty group to be used as the partition header
  ---------------------------------------------------------------------
*/
struct group *
CreateHeaderGroup()
{
  struct group *temp;

  temp = (struct group *)calloc(1,sizeof(struct group));

  temp->label = -1;
  temp->size = -1;
  temp->totlinks = -1;
  temp->inlinks = -1;
  temp->outlinks = -1;
  temp->totlinksW = -1.0;
  temp->inlinksW = -1.0;
  temp->outlinksW = -1.0;

  temp->coorX = -1.0;
  temp->coorY = -1.0;
  temp->coorZ = -1.0;

  temp->nodeList = NULL;
  temp->next = NULL;

  temp->offspr = NULL;

  return temp;
}

/*
  ---------------------------------------------------------------------
  Create a group at the end of the list starting at part
  ---------------------------------------------------------------------
*/
struct group *
CreateGroup(struct group *part, int label)
{
  while (part->next != NULL)
    part = part->next;

  part->next = (struct group *)calloc(1, sizeof(struct group));
  (part->next)->label = label;
  (part->next)->size = 0;
  (part->next)->totlinks = 0;
  (part->next)->inlinks = 0;
  (part->next)->outlinks = 0;
  (part->next)->totlinksW = 0.0;
  (part->next)->inlinksW = 0.0;
  (part->next)->outlinksW = 0.0;

  (part->next)->coorX = -1.0;
  (part->next)->coorY = -1.0;
  (part->next)->coorZ = -1.0;

  (part->next)->nodeList =
    (struct node_lis *)calloc(1, sizeof(struct node_lis));
  (part->next)->next = NULL;
  (part->next)->offspr = NULL;

  ((part->next)->nodeList)->node = -1;
  ((part->next)->nodeList)->status = -1;
  ((part->next)->nodeList)->next = NULL;
  ((part->next)->nodeList)->ref = NULL;
  ((part->next)->nodeList)->btw = -1.0;
  ((part->next)->nodeList)->weight = -1.0;

  return part->next;
}


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Partition creation
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/
/**
Read a partition from a file, the file should contain a line by
module, with the names of the nodes separated by tabulations.
 **/
struct group*
FReadPartition(FILE *inF){
  char label[MAX_LABEL_LENGTH];
  char sep[2];
  int noReadItems = 0;

  struct group *g = NULL;
  struct group *part = NULL;
  int npart = 0, nfields = 0;

  part = CreateHeaderGroup();
  g = CreateGroup(part, npart);

  while (!feof(inF)){
	nfields=fscanf(inF,"%[^\t\n]%[\t\n]",&label,&sep);
	if (nfields) {
	  AddNodeToGroupSoft(g, label);
	  if(sep[0]=='\n'){
		npart++;
		g = CreateGroup(part, npart);
	  }
	}
  }
  return(part);
}

/*
  ---------------------------------------------------------------------
  Build a partition from a file. The file should contain
  ---------------------------------------------------------------------
*/
struct group *
FCreatePartition(FILE *inF)
{
  char label[MAX_LABEL_LENGTH];
  char *separator = "///";
  struct group *g = NULL;
  struct group *part = NULL;
  int npart = 0;
  int noReadItems;

  /* Create the header of the partition */
  part = CreateHeaderGroup();

  /* Read the file */
  while (!feof(inF)) {
    g = CreateGroup(part, npart);
    npart++;
    noReadItems = fscanf(inF, "%s\n", &label[0]);
	if (noReadItems != 1)
	  printf ("Failed to read input, incorrect field number (%d != 1)\n",noReadItems);

    while (strcmp(label, separator) != 0) {
      AddNodeToGroupSoft(g, label);
      noReadItems = fscanf(inF, "%s\n", &label[0]);
	  if (noReadItems != 1)
		printf ("Failed to read input, incorrect field number (%d != 1)\n",noReadItems);

    }
  }

  return part;
}

/*
  ---------------------------------------------------------------------
  For an arbitrary network, create a partition with groups of nodes
  of a given size. The nodes are placed in groups in order.
  ---------------------------------------------------------------------
*/
struct group *
CreateEquiNPartition(struct node_gra *net, int gsize)
{
  int i,j;
  int ngroups;
  struct group *g = NULL;
  struct group *part = NULL;
  struct node_gra *p = NULL;

  /* Initialize stuff */
  p = net;
  part = CreateHeaderGroup();
  ngroups = CountNodes(net) / gsize;

  /* Create the groups and assign the nodes */
  for (i=0; i<=ngroups; i++) {
    g = CreateGroup(part,i);
    for (j=0; j<gsize; j++) {
      if ((p = p->next) != NULL) {
	AddNodeToGroup(g, p);
      }
    }
  }

  /* Done */
  return CompressPart(part);
}

/*
  ---------------------------------------------------------------------
  Create a partition with a given number of groups of a given
  size. The partition is not mapped to any network, and labels are
  set to integer numbers from 1 to N.
  ---------------------------------------------------------------------
*/
struct group *
CreateEquiNPartitionSoft(int ngroups, int gsize)
{
  int i, j;
  struct group *g = NULL;
  struct group *part = NULL;
  char label[MAX_LABEL_LENGTH];

  /* Initialize */
  part = CreateHeaderGroup();

  /* Create the groups and assign the nodes */
  for (i=0; i<ngroups; i++) {
    g = CreateGroup(part, i);
    for (j=0; j<gsize; j++) {
      sprintf(&label[0], "%d", (j + i*gsize + 1));
      AddNodeToGroupSoft(g, label);
    }
  }

  /* Done */
  return part;
}

/*
  ---------------------------------------------------------------------
  Given a network, create a partition from the inGroup attribute of
  its nodes
  ---------------------------------------------------------------------
*/
struct group *
CreatePartitionFromInGroup(struct node_gra *net)
{
  int i;
  int maxNGroup = 0;
  struct group **glist;
  struct group *part = NULL;
  struct node_gra *p = NULL;

  /* Determine the largest group number */
  p = net;
  while ((p = p->next) != NULL)
    if (p->inGroup > maxNGroup)
      maxNGroup = p->inGroup;

  /* Allocate enough memory for a list of pointers to groups (for
     faster access to groups), and initialize the list */
  glist = (struct group **) calloc(maxNGroup + 1, sizeof(struct group *));
  for(i=0; i<=maxNGroup; i++)
    glist[i] = NULL;

  /* Create the partition */
  part = CreateHeaderGroup();

  /* Assign the nodes creating new groups when necessary */
  p = net;
  while ((p = p->next) != NULL) {
    if (glist[p->inGroup] == NULL) {
      glist[p->inGroup] = CreateGroup(part, p->inGroup);
    }
    AddNodeToGroup(glist[p->inGroup], p);
  }

  /* Free memory and return */
  free(glist);
  return part;
}


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Partition removal
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Free the memory allocated to a partition
  ---------------------------------------------------------------------
*/
void
RemovePartition(struct group *part)
{
  if (part->next != NULL)
    RemovePartition(part->next);

  if (part->nodeList != NULL)
    FreeAdjacencyList(part->nodeList);

  if (part->offspr != NULL)
    RemovePartition(part->offspr);

  free(part);
}


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Node-group functions
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Add a node to a group, pointing all pointers and updating group and
  node properties.
  ---------------------------------------------------------------------
*/
struct node_lis *
AddNodeToGroup(struct group *g, struct node_gra *node)
{
  struct node_lis *p = g->nodeList;
  int totlink, inlink;
  double totweight, inweight;

  /* Go to the end of the list of nodes in the group */
  while (p->next != NULL)
    p = p->next;

  /* Create the node_lis and point it to the node */
  p->next = (struct node_lis *)calloc(1, sizeof(struct node_lis));
  (p->next)->node = node->num;
  (p->next)->nodeLabel = (char *) calloc(MAX_LABEL_LENGTH, sizeof(char));
  strcpy((p->next)->nodeLabel, node->label);
  (p->next)->status = 0;
  (p->next)->next = NULL;
  (p->next)->ref = node;

  (p->next)->btw=0.0;  /* initalising rush variable to 0.0; */
  (p->next)->weight=0.0;  /* initalising rush variable to 0.0; */

  /* Update the properties of the group */
  g->size++;
  totlink = NodeDegree(node);
  inlink = NLinksToGroup(node, g);
  totweight = NodeStrength(node);
  inweight = StrengthToGroup(node, g);
  g->totlinks += totlink - inlink;
  g->inlinks += inlink;
  g->outlinks = g->totlinks - g->inlinks;
  g->totlinksW += totweight - inweight;
  g->inlinksW += inweight;
  g->outlinksW = g->totlinksW - g->inlinksW;

  /* Update the properties of the node */
  node->inGroup = g->label;

  /* Done */
  return p->next;
}

/**
@brief Add a node to a graph without fully updating the group properties.
@author DB Stouffer

Like AddNodeToGroup() but the group properties totlinks, inlinks,
outlinks and their weighted counterpart are not updated.
**/
struct node_lis *
AddNodeToGroupFast(struct group *g, struct node_gra *node)
{
  struct node_lis *p = g->nodeList;

  /* Go to the end of the list of nodes in the group */
  while (p->next != NULL)
    p = p->next;

  /* Create the node_lis and point it to the node */
  p->next = (struct node_lis *)calloc(1, sizeof(struct node_lis));
  (p->next)->node = node->num;
  (p->next)->nodeLabel = (char *) calloc(MAX_LABEL_LENGTH, sizeof(char));
  strcpy((p->next)->nodeLabel, node->label);
  (p->next)->status = 0;
  (p->next)->next = NULL;
  (p->next)->ref = node;

  (p->next)->btw=0.0;  /* initalising rush variable to 0.0; */
  (p->next)->weight=0.0;  /* initalising rush variable to 0.0; */

  /* Partially update the properties of the group */
  g->size++;

  /* Update the properties of the node */
  node->inGroup = g->label;

  /* Done */
  return p->next;
}

/*
  ---------------------------------------------------------------------
  Softly add a node to a group. Pointers do NOT point anywhere and
  there is NO updating group or node properties.
  ---------------------------------------------------------------------
*/
struct node_lis *
AddNodeToGroupSoft(struct group *g, char *label)
{
  struct node_lis *p = g->nodeList;

  /* Go to the end of the list of nodes in the group */
  while (p->next != NULL)
    p = p->next;

  /* Create the node_lis */
  p->next = (struct node_lis *)calloc(1, sizeof(struct node_lis));
  (p->next)->node = -1;
  (p->next)->nodeLabel = (char *) calloc(MAX_LABEL_LENGTH, sizeof(char));
  strcpy((p->next)->nodeLabel, label);
  (p->next)->status = 0;
  (p->next)->next = NULL;
  (p->next)->ref = NULL;

  (p->next)->btw=0.0;  /* initalising rush variable to 0.0; */
  (p->next)->weight=0.0;  /* initalising rush variable to 0.0; */

  /* Update the properties of the group */
  g->size++;

  /* Done */
  return p->next;
}

/*
   ---------------------------------------------------------------------
   Remove a node from a group. Returns 1 if the node has been
   successfully removed and 0 if the node is not found in the group.
   ---------------------------------------------------------------------
*/
int
RemoveNodeFromGroup(struct group *g, struct node_gra *node)
{
  struct node_lis *p = g->nodeList;
  struct node_lis *temp;
  int totlink,inlink;
  double totweight,inweight;

  /* Find the node */
  while ((p->next != NULL) && ((p->next)->ref != node))
    p = p->next;

  if (p->next == NULL)
    return 0;
  else{
    temp = p->next;
    p->next = (p->next)->next;
    free(temp->nodeLabel);
    free(temp);

    /* Update the properties of the group */
    g->size--;
    totlink = NodeDegree(node);
    inlink = NLinksToGroup(node, g);
    totweight = NodeStrength(node);
    inweight = StrengthToGroup(node, g);
    g->totlinks -= totlink - inlink;
    g->inlinks -= inlink;
    g->outlinks = g->totlinks - g->inlinks;
    g->totlinksW -= totweight - inweight;
    g->inlinksW -= inweight;
    g->outlinksW = g->totlinksW - g->inlinksW;

    /* Done */
    return 1;
  }
}

/**
@brief Remove a node from a group without fully updating the group properties.
@author DB Stouffer

Returns 1 if the node has been successfully removed and 0 if the node
is not found in the group.

Like RemoveNodeFromGroup() but the group properties totlinks, inlinks,
outlinks and their weighted counterpart are not updated.
*/
int
RemoveNodeFromGroupFast(struct group *g, struct node_gra *node)
{
  struct node_lis *p = g->nodeList;
  struct node_lis *temp;

  /* Find the node */
  while ((p->next != NULL) && ((p->next)->ref != node))
    p = p->next;

  if (p->next == NULL)
    return 0;
  else{
    temp = p->next;
    p->next = (p->next)->next;
    free(temp->nodeLabel);
    free(temp);

    /* Partially update the properties of the group */
    g->size--;

    /* Done */
    return 1;
  }
}

/*
  ---------------------------------------------------------------------
  Move a node from old group to new group. Returns 1 if successful,
  and 0 if node is not in the old group.
  ---------------------------------------------------------------------
*/
int
MoveNode(struct node_gra *node, struct group *old, struct group *new)
{
  if (RemoveNodeFromGroup(old,node) == 0)
    return 0;

  AddNodeToGroup(new, node);
  return 1;
}

/**
@brief Move a node from a group without fully updating the group properties.
@author DB Stouffer

Returns 1 if the node has been successfully removed and 0 if the node
is not in the old group.

Like MoveNode() but the group properties totlinks, inlinks,
outlinks and their weighted counterpart are not updated.
*/
int
MoveNodeFast(struct node_gra *node, struct group *old, struct group *new)
{
  if (RemoveNodeFromGroupFast(old,node) == 0)
    return 0;

  AddNodeToGroupFast(new, node);
  return 1;
}

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Group and partition operations
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Removes all empty groups from a partition
  ---------------------------------------------------------------------
*/
struct group *
CompressPart(struct group *part)
{
  struct group *g1;
  struct group *gtemp;

  g1 = part;
  while (g1->next != NULL) {
    if (g1->next->size == 0) {
      FreeNodeLis(g1->next->nodeList);
      gtemp = g1->next;
      g1->next = g1->next->next;
      free(gtemp);
      gtemp = NULL;
    }

    else{
      g1 = g1->next;
    }
  }

  return part;
}

/*
  ---------------------------------------------------------------------
  Given a partition and a group label, return the group with the label
  ---------------------------------------------------------------------
*/
struct group *
GetGroup(struct group *part, int label)
{
  struct group *g=part;

  while ((g = g->next) != NULL)
    if (g->label == label)
      return g;

  return NULL;
}

/*
  ---------------------------------------------------------------------
  Count the number of groups (after a given group, usually the header
  of the partition)
  ---------------------------------------------------------------------
*/
int
NGroups(struct group *part)
{
  int ngroup = 0;

  while ((part = part->next) != NULL)
    ngroup++;

  return ngroup;
}

/*
  ---------------------------------------------------------------------
  Count the number of non-empty groups (after a given group, usually
  the header of the partition)
  ---------------------------------------------------------------------
*/
int
NNonEmptyGroups(struct group *part)
{
  int ngroup = 0;

  while ((part = part->next) != NULL)
    if (part->size > 0)
      ngroup++;

  return ngroup;
}

/*
  ---------------------------------------------------------------------
  Count the number of nodes in all groups (after a given group,
  usually the header of the partition)
  ---------------------------------------------------------------------
*/
int
PartitionSize(struct group *part)
{
  int nnod = 0;

  while ((part = part->next) != NULL)
    nnod += part->size;

  return nnod;
}

/*
  ---------------------------------------------------------------------
  Removes all links to and from the nodes in a group
  ---------------------------------------------------------------------
*/
void
RemoveWithinGroupLinks(struct group *g, int symmetric_sw)
{
  struct node_lis *nod;
  struct node_lis *nei;

  nod = g->nodeList;
  while ((nod = nod->next) != NULL) {
    nei = nod->ref->neig;
    while ((nei = nei->next) != NULL)
      RemoveLink(nod->ref, nei->ref, symmetric_sw);
  }

  return;
}

/*
  ---------------------------------------------------------------------
  Removes all links between groups
  ---------------------------------------------------------------------
*/
void
RemoveBetweenGroupLinks(struct group *part, int symmetric_sw)
{
  struct group *g = part;
  struct node_lis *nod;
  struct node_lis *nei;

  while ((g = g->next) != NULL) {
    nod = g->nodeList;
    while ((nod = nod->next) != NULL) {
      nei = nod->ref->neig;
      while ((nei = nei->next) != NULL) {
	if (nod->ref->inGroup != nei->ref->inGroup)
	  RemoveLink(nod->ref,nei->ref, symmetric_sw);
      }
    }
  }

  return;
}

/*
  ---------------------------------------------------------------------
  Given a partition, calculate the 'pattern of connections' between
  groups. 'Pattern of connections' type_sw can be: 'n' raw number of
  connections; 'f': fraction of all connections; 'p' probability of
  connection; 'e' fraction of connections minus expected fraction of
  connections (as in the modularity); 'z' z-score. CAUTION: MAKE SURE
  ALL CALCULATIONS ARE RIGHT BEFORE USING!!!!!!!!!!
  ---------------------------------------------------------------------
*/
void
BlockModel(struct group *part, char type_sw, int list_sw)
{
  struct group *g1, *g2;
  int links = 0;
  int nodes = 0;
  double prob;
  double bij=0, av, sig;

  /* Count the total number of nodes and links, and the average
     linking probability */
  g1 = part;
  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
      nodes += g1->size;
      links += g1->inlinks;
      g2 = g1;
      while ((g2 = g2->next) != NULL) {
	if (g2->size > 0) {
	  links += NG2GLinks(g1, g2);
	}
      }
    }
  }
  prob = 2.0 * (double)links / (double)(nodes * (nodes - 1));

  /* Calculate and print the blockmodel */
  g1 = part;
  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
      g2 = part;
      while ((g2 = g2->next) != NULL) {
	if (g2->size > 0) {

	  /* Calculate the matrix element */
	  switch (type_sw) {
	  case 'n':
	    if (g1 == g2)
	      bij = (double)g1->inlinks;
	    else
	      bij = (double)NG2GLinks(g1, g2);
	    break;
	  case 'f':
	    if (g1 == g2)
	      bij = (double)g1->inlinks / (double)links;
	    else
	      bij = (double)NG2GLinks(g1, g2) / (double)links;
	    break;
	  case 'p':
	    if (g1 == g2)
	      bij = (double)NG2GLinks(g1, g2) /
		(double)(g1->size * (g1->size - 1));
	    else
	      bij = (double)NG2GLinks(g1, g2) /
		(double)(g1->size * g2->size);
	    break;
	  case 'e':
	    if (g1 == g2)
	      bij = (double)g1->inlinks / (double)links -
		(double)((g1->inlinks + g1->totlinks) *
			 (g2->inlinks + g2->totlinks)) /
		((double)(4 * links * links));
	    else
	      bij = (double)NG2GLinks(g1, g2) / (double)links -
		(double)((g1->inlinks + g1->totlinks) *
			 (g2->inlinks + g2->totlinks)) /
		((double)(4 * links * links));
	    break;
	  case 'z':
	    if (g1 == g2) {
	      av = prob * (double)(g1->size * (g1->size - 1) / 2.0);
	      sig = sqrt((double)(g1->size * (g1->size - 1) / 2.0) *
			 prob * (1.0 - prob));
	      bij = ((double)g1->inlinks - av) / sig;
	    }
	    else {
	      av = prob * (double)(g1->size * g2->size);
	      sig = sqrt((double)(g1->size * g2->size) * prob * (1.0 - prob));
	      bij = ((double)NG2GLinks(g1, g2) - av) / sig;
	    }
	    break;
	  }
	  printf("%d %d %lf\n", g1->label+1, g2->label+1, bij);
	}
      }
    }
  }
}

/*
  ---------------------------------------------------------------------
  Count the number of links from a node to a given group
  ---------------------------------------------------------------------
*/
int
NLinksToGroup(struct node_gra* node, struct group *g)
{
  struct node_lis *nei = node->neig;
  int inlink = 0;

  while ((nei = nei->next) != NULL)
    if ((nei->ref)->inGroup == g->label)
      inlink++;

  return inlink;
}

/*
  ---------------------------------------------------------------------
  Count the number of links *with a certain weight* from a node to a
  given group
  ---------------------------------------------------------------------
*/
int
NWeightLinksToGroup(struct node_gra* node, struct group *g, double w)
{
  struct node_lis *nei = node->neig;
  int inlink = 0;

  while ((nei = nei->next) != NULL)
    if (nei->weight == w && (nei->ref)->inGroup == g->label)
      inlink++;

  return inlink;
}

/*
  ---------------------------------------------------------------------
  Count the number of links from a node to a given group, based on the
  label of the group only
  ---------------------------------------------------------------------
*/
int
NLinksToGroupByNum(struct node_gra* node, int gLabel)
{
  struct node_lis *nei = node->neig;
  int inlink = 0;

  while ((nei = nei->next) != NULL)
    if ((nei->ref)->inGroup == gLabel)
      inlink++;

  return inlink;
}




/*
  ---------------------------------------------------------------------
  Weight of links from a node to a given group
  ---------------------------------------------------------------------
*/
double
StrengthToGroup(struct node_gra* node, struct group *g)
{
  struct node_lis *nei = node->neig;
  double inlink = 0.0;

  while ((nei = nei->next) != NULL)
    if ((nei->ref)->inGroup == g->label)
      inlink += nei->weight;

  return inlink;
}

/**
Weight of links from a node to a given group, based on the label of
the group only.
**/
double
StrengthToGroupByNum(struct node_gra* node, int gLabel)
{
  struct node_lis *nei = node->neig;
  double inlink = 0.0;

  while ((nei = nei->next) != NULL)
    if ((nei->ref)->inGroup == gLabel)
      inlink += nei->weight;

  return inlink;
}

/*
  ---------------------------------------------------------------------
  Count the number of links between a pair of groups
  ---------------------------------------------------------------------
*/
int
NG2GLinks(struct group *g1, struct group *g2)
{
  struct node_lis *p = g1->nodeList;
  int nlink = 0;

  while ((p = p->next) != NULL)
    nlink += NLinksToGroup(p->ref, g2);

  if (g1 == g2)
    return nlink / 2;
  else
    return nlink;
}

/*
  ---------------------------------------------------------------------
  Count the number of links with a certain weight between a pair of
  groups
  ---------------------------------------------------------------------
*/
int
NWeightG2GLinks(struct group *g1, struct group *g2, double w)
{
  struct node_lis *p = g1->nodeList;
  int nlink = 0;

  while ((p = p->next) != NULL)
    nlink += NWeightLinksToGroup(p->ref, g2, w);

  if (g1 == g2)
    return nlink / 2;
  else
    return nlink;
}

/*
  ---------------------------------------------------------------------
  Total weight of links between a pair of groups
  ---------------------------------------------------------------------
*/
double
NG2GLinksWeight(struct group *g1, struct group *g2)
{
  struct node_lis *p = g1->nodeList;
  double nlink = 0.0;

  while ((p = p->next) != NULL)
    nlink += StrengthToGroup(p->ref, g2);

  if (g1 == g2)
    return nlink / 2.0;
  else
    return nlink;
}

/*
  ---------------------------------------------------------------------
  Move all nodes in group g1 to group g2
  ---------------------------------------------------------------------
*/
void
MergeGroups(struct group *g1, struct group *g2)
{
  struct node_lis *p = g1->nodeList;

  while (p->next != NULL)
    MoveNode((p->next)->ref, g1, g2);
  return;
}

/* a la DB Stouffer */
void
MergeGroupsFast(struct group *g1, struct group *g2)
{
  struct node_lis *p = g1->nodeList;

  while (p->next != NULL)
    MoveNodeFast((p->next)->ref, g1, g2);
  return;
}

/*
  ---------------------------------------------------------------------
  Creates a copy of a group and puts it at the end of the list of
  groups hanging from copy_root
  ---------------------------------------------------------------------
*/
struct group *
CopyGroup(struct group *copy_root, struct group *g)
{
  struct group *copy;
  struct node_lis *node;

  copy = CreateGroup(copy_root, g->label);

  /* Copy the trivial fields */
  copy->coorX = g->coorX;
  copy->coorY = g->coorY;
  copy->coorZ = g->coorZ;

  /* Copy the node list */
  node = g->nodeList;
  while (node->next != NULL) {
    node = node->next;
    AddNodeToGroup(copy, node->ref);
  }

  /* Some recursivity is in order here :( */
  /* AND NEEDS TO BE TESTED BEFORE USING!!!!!!!!!! */
  if ( g->offspr != NULL )
    copy->offspr = CopyPartition(g->offspr);

  return copy;
}

/*
  ---------------------------------------------------------------------
  Creates a copy of a whole partition. CAUTION!!!! The network MUST be
  mapped to the copy or the original partition if they are to be used
  after the copy has been created!!!!
  ---------------------------------------------------------------------
*/
struct group *
CopyPartition(struct group *original)
{
  struct group *copy_root = NULL;
  struct group *copy_p = NULL;
  struct group *part = original;

  copy_p = copy_root = CreateHeaderGroup();

  while ((part = part->next) != NULL) {
    CopyGroup(copy_p, part);
    copy_p = copy_p->next;
  }

  return copy_root;
}

/*
  ---------------------------------------------------------------------
  Creates a network that contains only the nodes in one of the groups
  of a partition
  ---------------------------------------------------------------------
*/
struct node_gra *
BuildNetFromGroup(struct group *group)
{
  struct node_gra *net = NULL;
  struct node_gra *last = NULL;
  struct node_gra *new = NULL;
  struct node_lis *p = (group->nodeList);

  last = net = CreateHeaderGraph();

  while ((p = p->next) != NULL) {
    new = CreateNodeGraph(last, p->nodeLabel);
    new->state = p->ref->state;
    new->coorX = p->ref->coorX;
    new->coorY = p->ref->coorY;
    new->coorZ = p->ref->coorZ;
    new->ivar1 = p->ref->ivar1;
    new->dvar1 = p->ref->dvar1;
    CopyAdjacencyList(p->ref, new);
    last = new;
  }

  CleanAdjacencies(net);
  RewireAdjacencyByLabel(net);
  RenumberNodes(net);

  return net;
}


/*
  ---------------------------------------------------------------------
  Creates a network that contains only the nodes in one of the groups
  of a partition AND its neighbors
  ---------------------------------------------------------------------
*/
struct node_gra *
BuildNetFromGroupNeig(struct group *group)
{
  struct node_gra *net = NULL;
  struct node_gra *last = NULL;
  struct node_gra *new = NULL;
  struct node_lis *p = group->nodeList;
  struct node_lis *nei = NULL;

  net = CreateHeaderGraph();
  last = net;

  /* Loop over nodes in the group */
  while ((p = p->next) != NULL){
    /* Add the node */
    new = CreateNodeGraph(last, p->ref->label);
    new->state = p->ref->state;
    new->coorX = p->ref->coorX;
    new->coorY = p->ref->coorY;
    new->coorZ = p->ref->coorZ;
    new->inGroup = p->ref->inGroup;
    last = new;
    CopyAdjacencyList(p->ref, new);

    /* Add the neighbors of the node */
    nei = p->ref->neig;
    while ((nei = nei->next) != NULL) {
      /* If the node does not exist and is in another module, add it
	 to the network */
      if ((IsThereNode(nei->ref->label, net) == 0) &&
	  (nei->ref->inGroup != p->ref->inGroup)) {
	new = CreateNodeGraph(last, nei->ref->label);
	new->state = nei->ref->state;
	new->coorX = nei->ref->coorX;
	new->coorY = nei->ref->coorY;
	new->coorZ = nei->ref->coorZ;
	new->inGroup = nei->ref->inGroup;
	last = new;
	CopyAdjacencyList(nei->ref, new);
      }
    }
  }

  CleanAdjacencies(net);
  RewireAdjacencyByNum(net);

  return net;
}

/*
  ---------------------------------------------------------------------
  Get some statistical properties (mean, std dev, min, and max) for
  the sizes of the modules in a partition.
  ---------------------------------------------------------------------
*/
void
GroupSizeStatistics(struct group *part,
		    double *theMean,
		    double *theStddev,
		    double *theMin,
		    double *theMax)
{
  int nGroups = NNonEmptyGroups(part);
  double *sizes = NULL;
  int count = 0;
  struct group *g = part;

  /* Allocate memory */
  sizes = allocate_d_vec(nGroups);

  /* Get group sizes */
  while ((g = g->next) != NULL) {
    if (g->size > 0) {
      sizes[count] = (double)g->size;
      count++;
    }
  }

  /* Get the statistical properties */
  (*theMean) = mean(sizes, nGroups);
  (*theStddev) = stddev(sizes, nGroups);
  (*theMin) = min(sizes, nGroups);
  (*theMax) = max(sizes, nGroups);

  /* Free memory and return */
  free_d_vec(sizes);
  return;
}


/*
  ---------------------------------------------------------------------
  Returns the first empty group in a partition. If no groups are empty,
  returns NULL.
  ---------------------------------------------------------------------
  */
struct group *
GetEmptyGroup(struct group *part)
{
  struct group *g;
  g = part;
  while ((g=g->next) != NULL) {
    if (g->size == 0) {
      return g;
    }
  }
  return NULL;
}



/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Network-partition operations
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Reset the inGroup attribute of all the nodes in a network
  ---------------------------------------------------------------------
*/
void
ResetNetGroup(struct node_gra *net)
{
  while ((net = net->next) != NULL)
    net->inGroup = -1;

  return;
}

/*
  ---------------------------------------------------------------------
  Given a partition and a network, set all the pointers and all the
  group and node attributes to the right values.
  ---------------------------------------------------------------------
*/
void
MapPartToNet(struct group *part, struct node_gra *net)
{
  struct group *g = NULL;
  struct node_lis *nod;
  int totlink, inlink;
  double totlinkW, inlinkW;
  void *nodeDict = NULL;
  struct node_gra *node = NULL;

  /* Reset the group of all nodes */
  ResetNetGroup(net);

  /* Create the nodeDict for fast access to nodes by label */
  nodeDict = MakeLabelDict(net);

  /* Go through the groups, reset attributes, and point pointers */
  g = part;
  while ((g = g->next) != NULL) {
    g->totlinks = 0;
    g->inlinks = 0;
    g->outlinks = 0;
    g->totlinksW = 0.0;
    g->inlinksW = 0.0;
    g->outlinksW = 0.0;

    nod = g->nodeList;
    while ((nod = nod->next) != NULL) {
      /* Get the node_gra by label */
      node = GetNodeDict(nod->nodeLabel, nodeDict);

      /* Update the properties of the group */
      nod->ref = node;
      nod->node = node->num;
      /* Update the properties of the node */
      node->inGroup = g->label;
    }
  }

  /* Count links in groups */
  g = part;
  while ((g = g->next) != NULL) {
    nod = g->nodeList;
    while ((nod = nod->next) != NULL) {
      totlink = NodeDegree(nod->ref);
      inlink = NLinksToGroup(nod->ref,g);
      totlinkW = NodeStrength(nod->ref);
      inlinkW = StrengthToGroup(nod->ref, g);
      g->totlinks += totlink;
      g->inlinks += inlink;
      g->totlinksW += totlinkW;
      g->inlinksW += inlinkW;
    }
    g->inlinks /= 2;
    g->totlinks -= g->inlinks;
    g->outlinks = g->totlinks - g->inlinks;
    g->inlinksW /= 2.0;
    g->totlinksW -= g->inlinksW;
    g->outlinksW = g->totlinksW - g->inlinksW;
  }

  /* Free memory allocated locally */
  FreeLabelDict(nodeDict);

  // Done
  return;
}

/* a la DB Stouffer */
void
MapPartToNetFast(struct group *part, struct node_gra *net)
{
  struct group *g = NULL;
  struct node_lis *nod;
  int totlink, inlink;
  double totlinkW, inlinkW;
  void *nodeDict = NULL;
  struct node_gra *node = NULL;

  /* Reset the group of all nodes */
  ResetNetGroup(net);

  /* Create the nodeDict for fast access to nodes by label */
  nodeDict = MakeLabelDict(net);

  /* Go through the groups, reset attributes, and point pointers */
  g = part;
  while ((g = g->next) != NULL) {
    g->totlinks = 0;
    g->inlinks = 0;
    g->outlinks = 0;
    g->totlinksW = 0.0;
    g->inlinksW = 0.0;
    g->outlinksW = 0.0;

    nod = g->nodeList;
    while ((nod = nod->next) != NULL) {
      /* Get the node_gra by label */
      node = GetNodeDict(nod->nodeLabel, nodeDict);

      /* Update the properties of the group */
      nod->ref = node;
      nod->node = node->num;
      /* Update the properties of the node */
      node->inGroup = g->label;
    }
  }

  /* Count links in groups */
/*  g = part;
  while ((g = g->next) != NULL) {
    nod = g->nodeList;
    while ((nod = nod->next) != NULL) {
      totlink = NodeDegree(nod->ref);
      inlink = NLinksToGroup(nod->ref,g);
      totlinkW = NodeStrength(nod->ref);
      inlinkW = StrengthToGroup(nod->ref, g);
      g->totlinks += totlink;
      g->inlinks += inlink;
      g->totlinksW += totlinkW;
      g->inlinksW += inlinkW;
    }
    g->inlinks /= 2;
    g->totlinks -= g->inlinks;
    g->outlinks = g->totlinks - g->inlinks;
    g->inlinksW /= 2.0;
    g->totlinksW -= g->inlinksW;
    g->outlinksW = g->totlinksW - g->inlinksW;
  }*/

  /* Free memory allocated locally */
  FreeLabelDict(nodeDict);

  // Done
  return;
}

/*
  ---------------------------------------------------------------------
  Same as MapPartToNet but ONLY the network (not the partition) is
  updated. This enables one to map a partition to networks that are a
  subset of the original network on which the partition is based.
  ---------------------------------------------------------------------
*/
void
MapPartToNetSoft(struct group *part, struct node_gra *net)
{
  struct group *g = NULL;
  struct node_lis *nod;
  void *nodeDict = NULL;
  struct node_tree *treeNode = NULL;
  struct node_gra *node = NULL;

  /* Reset the group of all nodes */
  ResetNetGroup(net);

  /* Create the nodeDict for fast access to nodes by label */
  nodeDict = MakeLabelDict(net);

  /* Go through the groups, reset attributes, and point pointers */
  g = part;
  while ((g = g->next) != NULL) {
    g->totlinks = 0;
    g->inlinks = 0;
    g->outlinks = 0;
    g->totlinksW = 0.0;
    g->inlinksW = 0.0;
    g->outlinksW = 0.0;

    nod = g->nodeList;
    while ((nod = nod->next) != NULL) {
      /* Get the node_gra by label */
      treeNode = CreateNodeTree();
      strcpy(treeNode->label, nod->nodeLabel);
      node = (*(struct node_tree **)tfind((void *)treeNode,
					  &nodeDict,
					  NodeTreeLabelCompare))->ref;
      FreeNodeTree(treeNode);
      /* Update the properties of the node */
      node->inGroup = g->label;
    }
  }

  return;
}

/*
  ---------------------------------------------------------------------
  Given an undirected network, create a partition in which each group
  corresponds to an isolated cluster of the network.
  ---------------------------------------------------------------------
*/
struct group *
ClustersPartition(struct node_gra *net)
{
  struct node_gra *p, *last;
  struct node_bfs *list,*lp,*lp2;
  int nnod=0;
  int anod=0;
  int size;
  int size_ant;
  struct group *part = NULL;
  struct group *thisgroup = NULL;
  int groupcoun = 0;
  int d;

  /* Initialize some variables */
  part = CreateHeaderGroup();
  size = 0;
  ResetNetGroup(net);
  nnod = CountNodes(net);
  list = CreateHeaderList();

  /* Start the search for clusters */
  last = net;
  do{
    /* Find the first unclassified node... */
    p = last->next;
    while (p->inGroup >= 0)
      p = p->next;
    last = p;

    /* Create a new group in the partition and add the node */
    thisgroup = CreateGroup(part, groupcoun++);
    AddNodeToGroup(thisgroup, p);

    /* ...and enqueue it */
    ResetNodesState(net);
    d = 0;
    Enqueue(p, list, list, &size, d);
    anod++;
    lp = list;

    /* Update the list with successive neighbors */
    do{
      size_ant = size;
      lp = RenewQueue(list, lp, &size, d);
      lp2 = lp;

      while (lp2->next != NULL) {
	lp2 = lp2->next;
	anod++;

	/* Add the node to the group */
	AddNodeToGroup(thisgroup, lp2->ref);
      }
      d++;
    } while(size != size_ant);

    ClearList(list, &size);
  } while(anod < nnod);

  /* Free memory */
  free(list);

  /* Done */
  return part;
}

/*
  ---------------------------------------------------------------------
  Remove all the links in a network that run between groups (as
  defined in the inGroup attribute of each node).
  ---------------------------------------------------------------------
*/
void
RemoveInterGroupLinks(struct node_gra *net)
{
  struct node_gra *p = net;
  struct node_lis *n;
  int thisgroup;

  while ((p = p->next) != NULL) {
    thisgroup = p->inGroup;
    n = p->neig;
    while ((n = n->next) != NULL) {
      if(n->ref->inGroup != thisgroup)
	RemoveLink(p, n->ref, 0);
    }
  }
}

/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Group and partition output
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Print the groups in a partition. If list_sw == 0, each group is
  printed in a row with some group info (nodes, links,
  etc.). Otherwise, each node is printed in a row and groups are
  separated by a line with '//'
  ---------------------------------------------------------------------
*/
void
FPrintPartition(FILE *outf, struct group *partition, int list_sw)
{
  struct group *g = partition;
  struct node_lis *p;

  while ((g = g->next) != NULL) {
    if(g->size > 0){

      /* Print group info for the verbose (non-list) mode */
      if (list_sw == 0) {
	fprintf(outf, "%d %d %d %d %d %lf %lf %lf ---",
		g->label+1, g->size, g->totlinks, g->inlinks, g->outlinks,
		g->totlinksW, g->inlinksW, g->outlinksW);
      }

      /* Print nodes */
      p = g->nodeList;
      while ((p = p->next) != NULL) {
	if (list_sw == 0)
	  fprintf(outf, " %s", p->nodeLabel);
	else
	  fprintf(outf, "%s\n", p->nodeLabel);
      }

      /* End group */
      if (list_sw == 0)
	fprintf(outf, "\n");
      else
	fprintf(outf, "///\n");
    }
  } /* End of loop over groups in the partition */

  return;
}

/*
  ---------------------------------------------------------------------
  Prints a partition to a file in Pajek format. The partition is
  printed from the inGroup attribute of the nodes of a network.
  ---------------------------------------------------------------------
*/
void
FPrintPajekPartitionFile(char *fname, struct node_gra *net)
{
  FILE *outF;

  outF = fopen(fname, "w");
  fprintf(outF, "*Vertices %d\n", CountNodes(net));

  while ((net = net->next) != NULL) {
    fprintf(outF, "   %d\n", net->inGroup + 1);
  }
  fclose(outF);

  return;
}


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Partition comparison
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Compute the mutual information between two partitions
  ---------------------------------------------------------------------
*/
double
MutualInformation(struct group *part1, struct group *part2, int label_sw)
{
  struct group *g1 = NULL, *g2 = NULL;
  struct node_lis *n1 = NULL, *n2 = NULL;
  int S = 0, S1 = 0, S2 = 0, S12 = 0;
  double H1 = 0.0, H2 = 0.0, H12 = 0.0;
  double I12 = 0.0;

  /*
    Count the number of nodes
  */
  S = PartitionSize(part1);
  if (S != PartitionSize(part2))
    fprintf(stderr, "WARNING : partitions have different size!\n");

  /*
    Compute the H1 and H2 entropies
  */
  g1 = part1;
  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
      S1 = g1->size;
      H1 += (double)S1 * log((double)S1 / (double)S);
    }
  }
  g2 = part2;
  while ((g2 = g2->next) != NULL) {
    if (g2->size > 0) {
      S2 = g2->size;
      H2 += (double)S2 * log((double)S2 / (double)S);
    }
  }

  /*
    Compute the join entropy H12
  */
  g1 = part1;
  while ((g1 = g1->next) != NULL) {
    if (g1->size > 0) {
      S1 = g1->size;

      g2 = part2;
      while ((g2 = g2->next) != NULL) {
	if (g2->size > 0) {

	  S2 = g2->size;

	  /*
	    Compute overlap
	  */
	  S12 = 0;
	  n1 = g1->nodeList;
	  while ((n1 = n1->next) != NULL) {
	    n2 = g2->nodeList;
	    while ((n2 = n2->next) != NULL) {
	      if ((label_sw == 0 && n1->node == n2->node ) || (label_sw == 1 && strcmp(n1->nodeLabel,n2->nodeLabel) == 0)) {
		S12++;
		break;
	      }
	    }
	  } /* Overlap size S12 has been calculated */
	  if (S12 > 0)
	    H12 += (double)S12 * log((double)(S12 * S) /
				 (double)(S1 * S2));
	}
      }
    }
  } /* End of loop over groups */

  /*
    Compute mutual information
  */
  I12 = -2.0 * H12 / (H1 + H2);

  return I12;
}

/*
  ---------------------------------------------------------------------
  Given a reference partition, calculate the fraction of nodes that
  are correctly classified in an actual partition. CAUTION: The
  function will not work if the largest group label is larger than the
  number of nodes in the network. But... why would that happen?!
  ---------------------------------------------------------------------
*/
double
CorrectlyClassified(struct group *refpart, struct group *actpart)
{
  int i;
  struct group *g1 = NULL;
  struct group *g2 = NULL;
  struct node_lis *p = NULL;
  int *group;
  int *score;
  int *max_score;
  int *map;
  int correct = 0;
  int nnod = PartitionSize(refpart);

  /* Allocate memory for all the arrays and initialize them*/
  group = allocate_i_vec(nnod);
  map = allocate_i_vec(nnod);
  score = allocate_i_vec(nnod);
  max_score = allocate_i_vec(nnod);
  for (i=0; i<nnod; i++) {
    group[i] = -1;
    map[i] = -1;
    score[i] = 0;
    max_score[i] = 0;
  }

  /* Find the reference group for each node */
  g1 = refpart;
  while ((g1 = g1->next) != NULL) {
    p = g1->nodeList;
    while ((p = p->next) != NULL) {
      group[p->node] = g1->label;
    }
  }

  /* Map the groups in the actual partition to the groups in the
     reference partition */
  g2 = actpart;
  while ((g2 = g2->next) != NULL) {
    p = g2->nodeList;
    while ((p = p->next) != NULL) {
      score[group[p->node]] += 1; /* update the scores */
    }

    /* pick the group with the highest score */
    g1 = refpart;
    while ((g1 = g1->next) != NULL) {
      if (score[g1->label] > max_score[g2->label]) {
	max_score[g2->label] = score[g1->label];
	map[g2->label] = g1->label;
      }
      score[g1->label] = 0;
    }
  }

  /* If two groups are mapped into the same group of the reference
     partition, get only the one that contains more nodes in that
     reference partition */
  g1 = actpart;
  while ((g1 = g1->next) != NULL) {
    g2 = g1;
    while ((g2 = g2->next) != NULL) {
      if (map[g1->label] == map[g2->label]){
	if (max_score[g1->label] > max_score[g2->label])
	  map[g2->label] = -1;
	else
	  map[g1->label] = -1;
      }
    }
  }

  /* Count correctly classified nodes */
  g2 = actpart;
  while ((g2 = g2->next) != NULL) {
    nnod += g2->size;
    if (map[g2->label] > -1){
      p = g2->nodeList;
      while ((p = p->next) != NULL) {
	if(group[p->node] == map[g2->label]){
	  correct++;
	}
      }
    }
  }

  /* Free memory allocated to the arrays */
  free_i_vec(group);
  free_i_vec(map);
  free_i_vec(score);
  free_i_vec(max_score);

  /* Done */
  return (double)correct / (double)nnod;
}


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Module indentification
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Calculate the modularity of a partition
  ---------------------------------------------------------------------
*/
double
Modularity(struct group *part)
{
  struct group *g = part;
  int links2 = 0;
  double modul = 0.0;

  /* Calculate the number of links times 2 */
  while ((g = g->next) != NULL)
    links2 += g->totlinks + g->inlinks;

  g = part;
  while ((g = g->next) != NULL)
    modul += (double)(2 * g->inlinks) / (double)links2 -
      ((double)(g->inlinks + g->totlinks) *
       (double)(g->inlinks + g->totlinks)) /
      ((double)links2 * (double)links2);

  return modul;
}

/*
  ---------------------------------------------------------------------
  Calculate the weighted modularity of a partition
  ---------------------------------------------------------------------
*/
double
ModularityWeight(struct group *part)
{
  struct group *g = part;
  double links2 = 0.0;
  double modul = 0.0;

  while ((g = g->next) != NULL)
    links2 += g->totlinksW + g->inlinksW;

  g = part;
  while (( g = g->next) != NULL)
    modul += (double)(2 * g->inlinksW) / (double)links2 -
      (g->inlinksW + g->totlinksW) *
      (g->inlinksW + g->totlinksW) /
      (links2 * links2);

  return modul;
}


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Roles
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Calculate the participation coefficient of a node
  ---------------------------------------------------------------------
*/
double
ParticipationCoefficient(struct node_gra *node)
{
  struct node_lis *nei = node->neig;
  int toGroup;
  double P = 0.0;
  int nlink = NodeDegree(node);

  // Go through the neighbors.
  if (nlink != 0) {
    while ((nei = nei->next) != NULL) {
      toGroup = NLinksToGroupByNum(node, nei->ref->inGroup);
      P += (double)toGroup / (double)(nlink * nlink);
    }
    P = 1.0 - P;
  }

  /* Done */
  return P;
}


/**

Compute the weighted participation coefficient of a node.
P_i = 1 - sum_{modules m}(strength_{im} / strength_i)**2

with:
- strength_i = sum of this node edges weigth.
- strength_{im} = sum of this node edges weigth.

**/
double
WeightedParticipationCoefficient(struct node_gra *node,  struct group *part)
{
  double toGroup;
  double P = 0.0;
  double strength = NodeStrength(node);
  double squared_st = strength*strength;
  //  printf ("%s: %f \n",node->label,strength);
  int nlink = NodeDegree(node);
  struct group *group=NULL;

  group = part;
  if (nlink != 0) {
	while ((group=group->next) != NULL){
      toGroup = StrengthToGroup(node, group);
      P += (toGroup*toGroup) / squared_st;
    }

    P = 1.0 - P;
  }
  //printf ("\n");

  return P;
}

/*
  ---------------------------------------------------------------------
  Calculate the within-module relative degree of a node. The network
  must be properly mapped to the partition under consideration.
  ---------------------------------------------------------------------
*/
double
WithinModuleRelativeDegree(struct node_gra *node, struct group *part)
{
  struct node_lis *p;
  double z;
  struct group *group=NULL;

  int inDegree;
  double kmean = 0.0, k2mean = 0.0, kstd;

  /* Find the group of the node */
  group = part;
  while (group->label != node->inGroup)
    group = group->next;

  /* Go through all the nodes in the group and calculate mean and
     standard deviation of the within-module degrees */
  p = group->nodeList;
  while ((p = p->next) != NULL) {
    inDegree = NLinksToGroup(p->ref, group);
    kmean += (double)inDegree;
    k2mean += (double)inDegree * (double)inDegree;
  }
  kmean /= (double)(group->size);
  k2mean /= (double)(group->size);
  kstd = sqrt(k2mean - kmean * kmean);

  /* Calculate the z-score */
  if (kstd == 0.0)
    z = 0.0;
  else
    z = ((double)NLinksToGroup(node, group) - kmean) / kstd;
  return z;
}




/**
Compute the the within-module relative strengh of a node.

Warning: The network must be properly mapped to the partition under
consideration.
**/
double
WithinModuleRelativeStrength(struct node_gra *node, struct group *part)
{
  struct node_lis *p;
  double z;
  struct group *group=NULL;

  double inStrength;
  double kmean = 0.0, k2mean = 0.0, kstd = 0.0;

  // Find the group of the node
  group = part;
  while (group->label != node->inGroup)
    group = group->next;

  // Go through all the nodes in the group and calculate mean and
  // standard deviation of the within-module strength
  p = group->nodeList;
  while ((p = p->next) != NULL) {
    inStrength = StrengthToGroup(p->ref, group);
    kmean += (double)inStrength;
    k2mean += (double)inStrength * (double)inStrength;
  }
  kmean /= (double)(group->size);
  k2mean /= (double)(group->size);
  kstd = sqrt(k2mean - kmean * kmean);

  // Calculate the z-score
  if (kstd == 0.0)
    z = 0.0;
  else
    z = ((double)StrengthToGroup(node, group) - kmean) / kstd;
  return z;
}


/*
  ---------------------------------------------------------------------
  Create a role partition according to the roles defined in Guimera &
  Amaral, Nature (2005). CAUTION: At the end of the process, the
  network is mapped onto the role partition.
  ---------------------------------------------------------------------
*/
struct group *
CatalogRoleIdent(struct node_gra *net, struct group *mod)
{
  struct group *roles = NULL;
  struct node_lis *p;
  int nroles = 7;
  int dest_group;
  int i;
  struct group *glist[7];
  struct group *g;
  double z, P;

  /* Create the groups */
  roles = CreateHeaderGroup();
  MapPartToNet(mod, net);
  glist[0] = CreateGroup(roles, 0);
  for (i=1; i<nroles; i++) {
    glist[i] = CreateGroup(glist[i-1],i);
  }

  /* Go through all the groups and assign roles to all the nodes */
  g = mod;
  while ((g = g->next) != NULL) {
    p = g->nodeList;
    while ((p = p->next) != NULL) {
      P = ParticipationCoefficient(p->ref);
      z = WithinModuleRelativeDegree(p->ref, g);
	  dest_group = GetRole(P,z);

      /* Add (softly) the node to the role group */
      AddNodeToGroupSoft(glist[dest_group], p->ref->label);
    } /* End of loop over nodes in this module */
  } /* End of loop over modules */

  /* Map the role partition onto the network and return*/
  MapPartToNet(roles, net);
  return roles;
}

/*
  ---------------------------------------------------------------------
  Create a role partition using the strength of the nodes rather than
  their degree.

  CAUTION: At the end of the process, the
  network is mapped onto the role partition.
  ---------------------------------------------------------------------
*/
struct group *
CatalogRoleIdentStrength(struct node_gra *net, struct group *mod)
{
  struct group *roles = NULL;
  struct node_lis *p;
  int nroles = 7;
  int dest_group;
  int i;
  struct group *glist[7];
  struct group *g;
  double z, P;

  /* Create the groups */
  roles = CreateHeaderGroup();
  MapPartToNet(mod, net);
  glist[0] = CreateGroup(roles, 0);
  for (i=1; i<nroles; i++) {
    glist[i] = CreateGroup(glist[i-1],i);
  }

  /* Go through all the groups and assign roles to all the nodes */
  g = mod;
  while ((g = g->next) != NULL) {
    p = g->nodeList;
    while ((p = p->next) != NULL) {
	  P = WeightedParticipationCoefficient(p->ref,mod);
	  z = WithinModuleRelativeStrength(p->ref, g);
	  dest_group = GetRole(P,z);

      /* Add (softly) the node to the role group */
      AddNodeToGroupSoft(glist[dest_group], p->ref->label);
    } /* End of loop over nodes in this module */
  } /* End of loop over modules */

  /* Map the role partition onto the network and return*/
  MapPartToNet(roles, net);
  return roles;
}




///////////////////////// SIMULATED ANNEALING //////////////////////////


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
      totallinks += NodeDegree(p);
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
	  nlink = NodeDegree(nodeList[target]);
	  
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
      totallinks += NodeDegree(p);
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
      totallinks += NodeDegree(p);
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
      nlink = NodeDegree(nlist[dice]);
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

/* RecursivePlotHier(FILE *nodes,FILE *links,FILE *sizes,struct group *part,int parentlabel,int *lastlabel,double scale) */
/* { */
/*   struct group *g = part; */

/*   while(g->next != NULL){ */
/*     g = g->next; */

/*     *lastlabel += 1; */

/*     fprintf(nodes,"%d \"\"      ellipse x_fact %g y_fact %g\n",*lastlabel,scale*sqrt(g->size),scale*sqrt(g->size)); */
/*     fprintf(links,"%d   %d   1\n",*lastlabel,parentlabel); */

/*     fprintf(sizes,"%d %d\n",*lastlabel,g->size); */

/*     if(g->offspr != NULL){ */
/*       RecursivePlotHier(nodes,links,sizes,g->offspr,*lastlabel,lastlabel,scale); */
/*     } */

/*   } */

/* } */

/* PlotPajekHierarchy(struct group *hier,double max_rad) */
/* { */
/*   int lastlabel = 1; */
/*   FILE *nodes,*links,*sizes; */
/*   int nnod = 0; */
/*   struct group *g = hier; */
/*   double scale; */

/*   while(g->next != NULL){ */
/*     g = g->next; */
/*     nnod += g->size; */
/*   } */
/*   scale = max_rad / sqrt(nnod); */

/*   nodes = fopen("vertices.net","w"); */
/*   links = fopen("links.net","w"); */
/*   sizes = fopen("sizes.dat","w"); */

/*   fprintf(sizes,"1 %d\n",nnod); */

/*   fprintf(links,"*Arcs\n*Edges\n"); */

/*   RecursivePlotHier(nodes,links,sizes,hier,1,&lastlabel,scale); */

/*   fclose(nodes); */
/*   fclose(links);   */
/*   fclose(sizes); */

/*   nodes = fopen("hier.net","w"); */
/*   fprintf(nodes,"*Vertices %d\n",lastlabel); */
/*   fprintf(nodes,"1 \"\"      ellipse x_fact %g y_fact %g\n",scale*sqrt(nnod),scale*sqrt(nnod)); */
/*   fclose(nodes); */
/* } */








/* struct group *ThermalNetworkSplitOdds(struct group *targ,double Ti, double Tf, gsl_rng *gen) */
/* { */
/*   struct group *glist[2]; */
/*   struct group *split = NULL; */
/*   struct node_gra *nlist[max_size]; */
/*   struct node_gra *net = NULL; */
/*   struct node_gra *p = NULL; */
/*   struct node_p; */
/*   int nnod = 0; */
/*   int i; */
/*   int des; */
/*   int target,oldg,newg; */
/*   int innew,inold,nlink; */
/*   int totallinks; */
/*   double dE; */
/*   double T, Ts = 0.95; */
/*   int A, B, C, D, Anew, Bnew, Cnew, Dnew; */

/*   for(i=0; i<max_size; i++){ */
/*     nlist[i] = NULL; */
/*   } */
/*   glist[0] = NULL; */
/*   glist[1] = NULL; */

/*   // Build a network from the nodes in the target group */
/*   net = BuildNetFromGroup(targ); */

/*   // Create the groups */
/*   split = CreateHeaderGroup(); */
/*   glist[0] = CreateGroup(split,0); */
/*   glist[1] = CreateGroup(split,1); */

/*   // Randomly assign the nodes to the groups */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     nlist[nnod] = p; */
/*     totallinks += NodeDegree(p); */
/*     nnod++; */

/*     des = floor(gsl_rng_uniform(gen)*2.0); */
/*     AddNodeToGroup(glist[des],p); */
/*   } */

/*   A = 2.0 * (glist[0]->size * glist[1]->size - glist[0]->outlinks); */
/*   B = 2.0 * glist[0]->outlinks; */
/*   C = glist[0]->size * (glist[0]->size - 1) + */
/*     glist[1]->size * (glist[1]->size - 1) - */
/*     2.0 * (glist[0]->inlinks + glist[1]->inlinks); */
/*   D = glist[0]->inlinks + glist[1]->inlinks; */

/*   // Do the SA to "optimize" the splitting */
/*   T = Ti; */
/*   while( T > Tf){ */

/*     for (i=0; i< nnod; i++){ */
/*       target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       if(oldg == 0) */
/* 	newg = 1; */
/*       else */
/* 	newg = 0; */

/*       // Calculate the change of energy */
/*       inold = NLinksToGroup(nlist[target],glist[oldg]); */
/*       innew = NLinksToGroup(nlist[target],glist[newg]); */
/*       nlink = NodeDegree(nlist[target]); */

/*       Anew = A - glist[newg]->size + innew + */
/* 	glist[oldg]->size - inold; */
/*       Bnew = B - innew + inold; */
/*       Cnew = C + glist[newg]->size - innew - */
/* 	glist[oldg]->size + inold; */
/*       Dnew = D + innew - inold;  */

/*       dE = Anew * Dnew / (Bnew * Cnew) - A * D / (B * C); */

/*       // Accept the change according to the Boltzman factor */
/*       if( (dE >= 0.0) || (gsl_rng_uniform(gen) < exp(dE/T)) ){ */
/* 	MoveNode(nlist[target],glist[oldg],glist[newg]); */
/*       } */
/*     } */

/*     T = T * Ts; */
/*   } */

/*   RemoveGraph(net); */

/*   return split; */
/* } */


/* // merge = 0 => No group merging */
/* struct group *SACommunityIdentOdds(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int merge, gsl_rng *gen) */
/* { */
/*   int i; */
/*   struct group *part = NULL; */
/*   struct group *split = NULL, *g = NULL; */
/*   struct group *glist[max_size]; */
/*   struct node_gra *nlist[max_size]; */
/*   struct node_gra *p; */
/*   struct node_lis *nod; */
/*   int target,empty; */
/*   int newg,oldg; */
/*   int nnod; */
/*   int totallinks = 0; */
/*   int innew,inold,nlink; */
/*   double energy = 0.0, dE; */
/*   double T; */
/*   int g1,g2; */
/*   double energyant = 0.0; */
/*   int count = 0, limit = 25; // to stop the search if the energy */
/*                           // does not change */
/*   int cicle1,cicle2; */
/*   int trans[maxim_int]; */
/*   int A, B, C, D, Anew, Bnew, Cnew, Dnew; */


/*   // Create the groups and assign each node to one group */
/*   nnod = CountNodes(net); */
/*   part = CreateHeaderGroup(); */
/*   p = net->next; */
/*   ResetNetGroup(net); // All nodes reset to group -1 */

/*   nlist[0] = p; */
/*   trans[p->num] = 0; */
/*   glist[0] = CreateGroup(part,0); */
/*   AddNodeToGroup(glist[0],p); */
/*   totallinks += NodeDegree(p); */

/*   for( i=1; i<nnod; i++ ) { */
/*     p = p->next; */

/*     nlist[i] = p; */
/*     trans[p->num] = i; */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*     AddNodeToGroup(glist[i],p); */
/*     totallinks += NodeDegree(p); */
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

/*   // Do the simulated annealing */
/*   T = Ti; */
/*   A = (nnod * (nnod -1) - totallinks) / 2; */
/*   B = totallinks / 2; */
/*   C = 0; */
/*   D = 0; */
/*   energy = 0.0; */

/*   while( T > Tf && count < limit){ */

/*     if (energy == energyant) */
/*       count++; */
/*     else{ */
/*       energyant = energy; */
/*       count = 0; */
/*     } */

/* /\*     PrintGroups(part); *\/ */
/* /\*     printf("%d %d %d %d\n",A,B,C,D); *\/ */
/* /\*     printf("%g %lf %lf %g\n",1.0/T, energy, Modularity(part), T); *\/ */
/*     printf("%g %lf %g\n",1.0/T, energy, T); */

/*     for ( i=0; i < cicle1; i++ ){ */

/*       /////////////////////////////// */
/*       // Propose an individual change */
/*       /////////////////////////////// */
/*       target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       }while(newg == oldg); */

/*       // Calculate the change of energy */
/*       inold = NLinksToGroup(nlist[target],glist[oldg]); */
/*       innew = NLinksToGroup(nlist[target],glist[newg]); */
/*       nlink = NodeDegree(nlist[target]); */

/*       Anew = A + (-glist[newg]->size + innew + */
/* 			glist[oldg]->size - 1 - inold); */
/*       Bnew = B - (innew - inold); */
/*       Cnew = C - (-glist[newg]->size + innew + */
/* 			glist[oldg]->size - 1 - inold); */
/*       Dnew = D + (innew - inold);  */

/*       if(B == 0 || C == 0) */
/* 	dE = 1.0; */
/*       else */
/* 	dE = (double)(Anew * Dnew) / (double)(Bnew * Cnew) - */
/* 	  (double)(A * D) / (double)(B * C); */

/* /\*       printf("%d %d %d %lf\n",innew,inold,nlink,exp(dE/T)); *\/ */

/*       // Accept the change according to Metroppolis */
/*       if( (dE >= 0.0) || ( gsl_rng_uniform(gen) < exp(dE/T) ) ){ */
/* 	MoveNode(nlist[target],glist[oldg],glist[newg]); */
/* 	energy = (double)(Anew * Dnew) / (double)(Bnew * Cnew); */
/* 	A = Anew; */
/* 	B = Bnew; */
/* 	C = Cnew; */
/* 	D = Dnew; */
/* /\* 	printf("  accepted %lf\n",energy); *\/ */
/*       } */
/*     } */

/*     if (merge == 1) { */

/*       for ( i=0; i < cicle2; i++ ){ */

/* 	////////////////////////////////////////////////////// */
/* 	// Propose a pair of merge/split collective changes // */
/* 	////////////////////////////////////////////////////// */

/* 	// Merge ///////////////////////////////////////////// */
/* 	target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/* 	g1 = nlist[target]->inGroup; */

/* 	if(glist[g1]->size < nnod){ */

/* 	  do{ */
/* 	    target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/* 	    g2 = nlist[target]->inGroup; */
/* 	  }while( g1 == g2 ); */

/* 	  // Calculate the change of energy */
/* 	  nlink = NG2GLinks(glist[g1],glist[g2]); */

/* 	  Anew = A - 2.0 * ( glist[g1]->totlinks * */
/* 			     glist[g2]->totlinks - nlink ); */
/* 	  Bnew = B - 2.0 * nlink; */
/* 	  Cnew = C + 2.0 * ( glist[g1]->totlinks * */
/* 			     glist[g2]->totlinks - nlink ); */
/* 	  Dnew = D + 2.0 * nlink; */

/* 	  if(B == 0 || C == 0) */
/* 	    dE = 1.0; */
/* 	  else */
/* 	    dE = Anew * Dnew / (Bnew * Cnew) - A * D / (B * C); */

/* 	  // Accept the change according to Metroppolis */
/* 	  if( (dE >= 0.0) || ( gsl_rng_uniform(gen) < exp(dE/T) ) ){ */
/* 	    MergeGroups(glist[g1],glist[g2]); */
/* 	    energy += dE; */
/* 	    A = Anew; */
/* 	    B = Bnew; */
/* 	    C = Cnew; */
/* 	    D = Dnew; */
/* 	  } */
/* 	} */

/* 	// Split ///////////////////////////////////////////// */
/* 	target = floor(gsl_rng_uniform(gen) * (double)nnod); // target node */
/* 	target = nlist[target]->inGroup;    // target group */

/* 	// Look for an empty group */
/* 	g = part; */
/* 	empty = -1; */
/* 	while((g->next != NULL) && (empty < 0)){ */
/* 	  g = g->next; */
/* 	  if (g->size == 0){ */
/* 	    empty = g->label; */
/* 	  } */
/* 	} */

/* 	if (empty >= 0 ){ // if there are no empty groups, do nothing */
/* 	  split = ThermalNetworkSplitOdds(glist[target],Ti,T,gen); */

/* 	  // Split the group */
/* 	  nod = (split->next)->nodeList; */
/* 	  while ( nod->next != NULL ){ */
/* 	    nod = nod->next; */
/* 	    MoveNode(nlist[trans[nod->node]], */
/* 		     glist[target], */
/* 		     glist[empty]); */
/* 	  } */
/* 	  RemovePartition(split); */
/* 	  split = NULL; */

/* 	  // Try to re-merge the two groups */
/* 	  // Calculate the change of energy */
/* 	  // Calculate the change of energy */
/* 	  nlink = NG2GLinks(glist[g1],glist[g2]); */

/* 	  Anew = A - 2.0 * ( glist[g1]->totlinks * */
/* 			     glist[g2]->totlinks - nlink ); */
/* 	  Bnew = B - 2.0 * nlink; */
/* 	  Cnew = C + 2.0 * ( glist[g1]->totlinks * */
/* 			     glist[g2]->totlinks - nlink ); */
/* 	  Dnew = D + 2.0 * nlink; */

/* 	  dE = Anew * Dnew / (Bnew * Cnew) - A * D / (B * C); */

/* 	  // Accept the change according to "inverse" Metroppolis. */
/* 	  // Inverse means that the algor is applied to the split */
/* 	  // and NOT to the merge! */
/* 	  if( (dE >= 0.0) && ( gsl_rng_uniform(gen) > exp(-dE/T) ) ){ */
/* 	    MergeGroups(glist[target],glist[empty]); */
/* 	  } */
/* 	  else{ */
/* 	    energy -= dE; */
/* 	    A = Anew; */
/* 	    B = Bnew; */
/* 	    C = Cnew; */
/* 	    D = Dnew; */
/* 	  } */

/* 	} */

/*       } */
/*     } */

/*     T = T * Ts; */
/*   } */

/*   printf("energy = %lf\n",energy); */
/*   printf("%d %d %d %d\n",A,B,C,D); */

/*   return CompressPart(part); */
/* } */


/* double BlockModularityE(struct group *part,int eij[],int nnod,double alpha) */
/* { */
/*   struct group *g = part; */
/*   int links = 0; */
/*   double role = 0.0; */
/*   struct group *glist[max_size]; */
/*   int i,j; */

/*   while ( g->next != NULL ){ */
/*     g = g->next; */
/*     glist[g->label] = g; */
/*     links += g->totlinks + g->inlinks; */
/*   } */

/*   for ( i=0; i<nnod; i++ ){ */

/*     if (glist[i]->size > 0) */
/*       for ( j=i; j<nnod; j++ ){ */

/* 	if (glist[j]->size > 0) */
/* 	  role += pow( */
/* 		      pow( eij[Map2Array(i,j,nnod)] / (double)links - */
/* 			   (glist[i]->totlinks + (glist[i]->inlinks)) * */
/* 			   (glist[j]->totlinks + (glist[j]->inlinks)) / */
/* 			   (double)(links * links) , 2.0 ), */
/* 		      alpha); */
/*       } */
/*   } */

/*   return role; */
/* } */

/* double BlockModularity(struct group *part,double alpha) */
/* { */
/*   struct group *g = part; */
/*   int links = 0; */
/*   double role = 0.0; */
/*   struct group *glist[max_size]; */
/*   int i,j; */
/*   int eij[max_size*(max_size+1)/2]; */
/*   int ngroup = 0; */

/*   // Map the groups to a list */
/*   while ( g->next != NULL ){ */
/*     g = g->next; */
/*     glist[g->label] = g; */
/*     links += g->totlinks + g->inlinks; */
/*     ngroup++; */
/*   } */

/*   // Calculate eij */
/*   for( i=0; i<ngroup; i++ ) { */
/*     eij[Map2Array(i,i,ngroup)] = 2 * glist[i]->inlinks; */
/*     for( j=i+1; j<ngroup; j++ ) { */
/*       eij[Map2Array(i,j,ngroup)] = 2 * */
/* 	NG2GLinks(glist[i],glist[j]); */
/*     } */
/*   } */

/*   // Calculate the "role energy" */
/*   for ( i=0; i<ngroup; i++ ){ */

/*     if (glist[i]->size > 0) */
/*       for ( j=i; j<ngroup; j++ ){ */

/* 	if (glist[j]->size > 0) */
/* 	  role += pow( */
/* 		      pow( eij[Map2Array(i,j,ngroup)] / (double)links - */
/* 			   (glist[i]->totlinks + (glist[i]->inlinks)) * */
/* 			   (glist[j]->totlinks + (glist[j]->inlinks)) / */
/* 			   (double)(links * links) , 2.0 ), */
/* 		      alpha); */
/*       } */
/*   } */

/*   return role; */
/* } */


/* struct group *SABlockModules(struct node_gra *net,double Ti,double Tf,double Ts,double fac,double alpha, gsl_rng *gen) */
/* { */
/*   int i,j,k; */
/*   struct group *part = NULL; */
/*   struct group *split = NULL, *g = NULL; */
/*   struct group *glist[max_size]; */
/*   struct node_gra *nlist[max_size]; */
/*   struct node_gra *p; */
/*   struct node_lis *nod; */
/*   int target,empty; */
/*   int newg,oldg; */
/*   int nnod; */
/*   int totallinks = 0; */
/*   double energy = 0.0, dE; */
/*   double T; */
/*   int g1,g2; */
/*   int eij[max_size*(max_size+1)/2]; */
/*   double eold,enew; */
/*   double energyant; */
/*   int count = 0, limit = 15; */

/*   // Create the groups and assign each node to one group */
/*   nnod = CountNodes(net); */
/*   part = CreateHeaderGroup(); */
/*   p = net->next; */

/*   nlist[0] = p; */
/*   glist[0] = CreateGroup(part,0); */
/*   AddNodeToGroup(glist[0],p); */
/*   totallinks += NodeDegree(p); */

/*   for( i=1; i<nnod; i++ ) { */
/*     p = p->next; */

/*     nlist[i] = p; */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*     AddNodeToGroup(glist[i],p); */
/*     totallinks += NodeDegree(p); */
/*   } */

/*   // Calculate the initial value of eij */
/*   for( i=0; i<nnod; i++ ) { */
/*     eij[Map2Array(i,i,nnod)] = 2 * glist[i]->inlinks; */
/*     for( j=i+1; j<nnod; j++ ) { */
/*       eij[Map2Array(i,j,nnod)] = 2 * */
/* 	NG2GLinks(glist[i],glist[j]); */
/*     } */
/*   } */

/*   // Do the simulated annealing */
/*   T = Ti; */
/*   energy = BlockModularityE(part,eij,nnod,alpha); */
/*   energyant = 0.0; */
/*   count = 0; */

/*   while( T > Tf && count < limit){ */

/*     if (fabs(energyant - energy) / fabs(energy) < 1.e-6) */
/*       count++; */
/*     else{ */
/*       energyant = energy; */
/*       count = 0; */
/*     } */

/*     printf("%g %lf %g\n",1.0/T, energy, T); */
/* /\*     printf("%g %lf %lf %g\n",1.0/T, energy, BlockModularityE(part,eij,nnod,alpha), T); *\/ */
/* /\*     printf("%g %lf %lf %g\n",1.0/T, energy, BlockModularity(part,alpha), T); *\/ */

/*     for ( k=0; k < (fac*nnod*nnod); k++ ){ */

/*       /////////////////////////////// */
/*       // Propose an individual change */
/*       /////////////////////////////// */
/*       target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       }while(newg == oldg); */

/*       // Calculate the change of energy */
/*       // Old energy */
/*       eold = enew = 0.0; */
/*       for( i=0; i<nnod; i++ ) { */
/* 	eold += pow( */
/* 		    pow( eij[Map2Array(oldg,i,nnod)] / */
/* 			 (double)totallinks - */
/* 			 (glist[oldg]->totlinks + */
/* 			  (glist[oldg]->inlinks)) * */
/* 			 (glist[i]->totlinks + (glist[i]->inlinks)) / */
/* 			 ((double)totallinks * (double)totallinks) , 2.0 ), */
/* 		    alpha); */
/* 	eold += pow( */
/* 		    pow( eij[Map2Array(newg,i,nnod)] / */
/* 			 (double)totallinks - */
/* 			 (glist[newg]->totlinks + */
/* 			  (glist[newg]->inlinks)) * */
/* 			 (glist[i]->totlinks + (glist[i]->inlinks)) / */
/* 			 ((double)totallinks * (double)totallinks) , 2.0 ), */
/* 		    alpha); */
/*       } */
/*       eold -= pow( */
/* 		  pow( eij[Map2Array(newg,oldg,nnod)] / */
/* 		       (double)totallinks - */
/* 		       (glist[newg]->totlinks + */
/* 			(glist[newg]->inlinks)) * */
/* 		       (glist[oldg]->totlinks + */
/* 			(glist[oldg]->inlinks)) / */
/* 		       ((double)totallinks * (double)totallinks) , 2.0 ), */
/* 		  alpha); */

/*       // Move the node and actualize the eij matrix */
/*       MoveNode(nlist[target],glist[oldg],glist[newg]); */
/*       nod = nlist[target]->neig; */
/*       while(nod->next != NULL){ */
/* 	nod = nod->next; */
/* 	eij[Map2Array(oldg,(nod->ref)->inGroup,nnod)] -= 2; */
/* 	eij[Map2Array(newg,(nod->ref)->inGroup,nnod)] += 2; */
/*       } */

/*       // New energy */
/*       for( i=0; i<nnod; i++ ) { */
/* 	enew += pow( */
/* 		    pow( eij[Map2Array(oldg,i,nnod)] / */
/* 			 (double)totallinks - */
/* 			 (glist[oldg]->totlinks + */
/* 			  (glist[oldg]->inlinks)) * */
/* 			 (glist[i]->totlinks + (glist[i]->inlinks)) / */
/* 			 ((double)totallinks * (double)totallinks) , 2.0 ), */
/* 		    alpha); */
/* 	enew += pow( */
/* 		    pow( eij[Map2Array(newg,i,nnod)] / */
/* 			 (double)totallinks - */
/* 			 (glist[newg]->totlinks + */
/* 			  (glist[newg]->inlinks)) * */
/* 			 (glist[i]->totlinks + (glist[i]->inlinks)) / */
/* 			 ((double)totallinks * (double)totallinks) , 2.0 ), */
/* 		    alpha); */
/*       } */
/*       enew -= pow( */
/* 		  pow( eij[Map2Array(newg,oldg,nnod)] / */
/* 		       (double)totallinks - */
/* 		       (glist[newg]->totlinks + */
/* 			(glist[newg]->inlinks)) * */
/* 		       (glist[oldg]->totlinks + */
/* 			(glist[oldg]->inlinks)) / */
/* 		       ((double)totallinks * (double)totallinks) , 2.0 ), */
/* 		  alpha); */

/*       // Accept the change according to Metroppolis */
/*       dE = enew - eold; */
/*       if( (dE >= 0.0) || ( gsl_rng_uniform(gen) < exp(dE/T) ) ){ */
/* 	energy += dE; */
/*       } */
/*       else{ */
/* 	// "Actualize back" the eij matrix and move the node back */
/* 	MoveNode(nlist[target],glist[newg],glist[oldg]); */
/* 	nod = nlist[target]->neig; */
/* 	while(nod->next != NULL){ */
/* 	  nod = nod->next; */
/* 	  eij[Map2Array(oldg,(nod->ref)->inGroup,nnod)] += 2; */
/* 	  eij[Map2Array(newg,(nod->ref)->inGroup,nnod)] -= 2; */
/* 	} */
/*       } */
/*     } */

/*     T = T * Ts; */
/*   } */

/*   printf("BlockModularity = %lf\n",BlockModularity(part,alpha)); */

/*   return part; */
/* } */

/* double BlockModularity3E(struct group *part,int eij[],int ngroup, int nlink, int nnod) */
/* { */
/*   struct group *g1, *g2; */
/*   double block = 0.0; */
/*   int nij, kij; */
/*   int **degen; */
/*   int Ng, max, dim; */
/*   int i, j; */
/*   int npair; */
/*   double p; */

/*   // Initialize the degeneration matrix */
/*   Ng = CountNonEmptyGroups(part); */
/*   max = SizeLargestGroup(part); */
/*   dim = max * max; */

/*   degen = allocate_i_mat(dim+1,dim+1); */

/*   for(i=0; i<=dim; i++){ */
/*     for(j=0; j<=dim; j++){ */
/*       degen[i][j] = 0; */
/*     } */
/*   } */

/*   // disregard links and pairs within the same block */
/*   npair = nnod * (nnod-1) / 2; */
/*   g1 = part; */
/*   while(g1->next != NULL){ */
/*     g1 = g1->next; */

/*     nlink -= g1->inlinks; */
/*     npair -= g1->size * (g1->size-1) / 2; */
/*   } */

/*   p = (double)nlink / (double)npair; */

/*   // Calculate the block modularity */
/*   g1 = part; */
/*   while(g1->next != NULL){ */
/*     g1 = g1->next; */

/*     if (g1->size > 0){ */

/*       // off-diagonal */
/*       g2 = g1; */
/*       while(g2->next != NULL){ */
/* 	g2 = g2->next; */

/* 	if (g2->size > 0){ */
/* 	  nij = g1->size * g2->size; */
/* 	  kij = eij[Map2Array(g1->label,g2->label,ngroup)]; */
/* 	  degen[nij][kij] += 1; */

/* 	  block += (double)nij * log( 1.0 - p ) + LogFactStir(nij) - */
/* 	    LogFactStir(nij-kij) - LogFactStir(kij); */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   // Add the link term to the modularity */
/*   if (nlink > 0) */
/*     block += (double)nlink * log(p/(1.0-p)); */

/*   // Add the degeneration term to the block modularity */
/*   block += LogFactStir(Ng*(Ng-1)/2); */

/*   for(i=0; i<dim; i++){ */
/*     for(j=0; j<dim; j++){ */
/*       if (degen[i][j] > 0){ */
/* 	block -= LogFactStir(degen[i][j]); */
/*       } */
/*     } */
/*   } */

/*   free_i_mat(degen,dim+1); */

/*   return block; */
/* } */


/* double BlockModularity3(struct group *part, struct node_gra *net) */
/* { */
/*   struct group *g1, *g2; */
/*   double block = 0.0; */
/*   int nij, kij; */
/*   int eij[max_size*(max_size+1)/2]; */
/*   int ngroup, nnod, nlink, npair; */
/*   double p; */
/*   int coun = 0; */
/*   int **degen; */
/*   int Ng, max, dim; */
/*   int i, j; */

/*   // Initialize the degeneration matrix */
/*   Ng = CountNonEmptyGroups(part); */
/*   max = SizeLargestGroup(part); */
/*   dim = max * max; */

/*   degen = allocate_i_mat(dim+1,dim+1); */

/*   for(i=0; i<dim; i++){ */
/*     for(j=0; j<dim; j++){ */
/*       degen[i][j] = 0; */
/*     } */
/*   } */

/*   ngroup = CountGroups(part); */
/*   nnod = CountNodes(net); */

/*   // disregard links and pairs within the same block */
/*   nlink = TotalNLinks(net); */
/*   npair = nnod * (nnod-1) / 2; */
/*   g1 = part; */
/*   while(g1->next != NULL){ */
/*     g1 = g1->next; */

/*     nlink -= g1->inlinks; */
/*     npair -= g1->size * (g1->size-1) / 2; */
/*   } */

/*   p = (double)nlink / (double)npair; */

/*   // Calculate eij */
/*   g1 = part; */
/*   while(g1->next != NULL){ */
/*     g1 = g1->next; */
/*     eij[Map2Array(g1->label,g1->label,ngroup)] = g1->inlinks; */

/*     g2 = g1; */
/*     while(g2->next != NULL){ */
/*       g2 = g2->next; */
/*       eij[Map2Array(g1->label,g2->label,ngroup)] = */
/* 	NG2GLinks(g1,g2); */
/*     } */
/*   } */

/*   // Calculate the block modularity */
/*   g1 = part; */
/*   while(g1->next != NULL){ */
/*     g1 = g1->next; */

/*     if (g1->size > 0){ */

/*       // off-diagonal */
/*       g2 = g1; */
/*       while(g2->next != NULL){ */
/* 	g2 = g2->next; */

/* 	if (g2->size > 0){ */
/* 	  nij = g1->size * g2->size; */
/* 	  kij = eij[Map2Array(g1->label,g2->label,ngroup)]; */
/* 	  degen[nij][kij] += 1; */

/* 	  block += (double)nij * log( 1.0 - p ) + LogFactStir(nij) - */
/* 	    LogFactStir(nij-kij) - LogFactStir(kij); */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   // Add the link term to the modularity */
/*   if (nlink > 0) */
/*     block += (double)nlink * log(p/(1.0-p)); */

/*   // Add the degeneration term to the block modularity */
/*   block += LogFactStir(Ng*(Ng-1)/2); */

/*   for(i=0; i<dim; i++){ */
/*     for(j=0; j<dim; j++){ */
/*       if (degen[i][j] > 0) */
/* 	block -= LogFactStir(degen[i][j]); */
/*     } */
/*   } */

/*   free_i_mat(degen,dim+1); */

/*   return block; */
/* } */


/* struct group *SABlockModules3(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int ngroup, gsl_rng *gen) */
/* { */
/*   int i,j,k; */
/*   struct group *part = NULL; */
/*   struct group *split = NULL, *g = NULL; */
/*   struct group *glist[max_size]; */
/*   struct node_gra *nlist[max_size]; */
/*   struct node_gra *p; */
/*   struct node_lis *nod; */
/*   int target,empty; */
/*   int newg,oldg; */
/*   int nnod; */
/*   int totallinks = 0; */
/*   double energy = 0.0, dE; */
/*   double T; */
/*   int g1,g2; */
/*   int eij[max_size*(max_size+1)/2]; */
/*   double eold,enew; */
/*   double energyant; */
/*   int count = 0, limit = 15; */
/*   int nlink; */

/*   nnod = CountNodes(net); */
/*   nlink = TotalNLinks(net); */

/*   // Create the groups  */
/*   part = CreateHeaderGroup(); */

/*   glist[0] = CreateGroup(part,0); */
/*   for( i=1; i<ngroup; i++ ) { */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*   } */

/*   // Assign nodes to groups at random */
/*   p = net; */
/*   for (i=0; i<nnod; i++) { */
/*     p = p->next; */
/*     nlist[i] = p; */
/*     totallinks += NodeDegree(p); */

/*     target = floor(gsl_rng_uniform(gen) * (double)ngroup); */
/*     AddNodeToGroup(glist[target],p); */
/*   } */

/*   // Calculate the initial value of eij */
/*   for( i=0; i<ngroup; i++ ) { */
/*     eij[Map2Array(i,i,ngroup)] = glist[i]->inlinks; */
/*     for( j=i+1; j<ngroup; j++ ) { */
/*       eij[Map2Array(i,j,ngroup)] = NG2GLinks(glist[i],glist[j]); */
/*     } */
/*   } */

/*   // Do the simulated annealing */
/*   T = Ti; */
/*   energy = BlockModularity3E(part,eij,ngroup,nlink,nnod); */
/*   energyant = 0.0; */
/*   count = 0; */

/*   while( T > Tf && count < limit){ */

/*     if (energyant == energy) */
/*       count++; */
/*     else{ */
/*       energyant = energy; */
/*       count = 0; */
/*     } */

/*     printf("%g %lf %g\n",1.0/T, energy, T); */
/* /\*     printf("%g %lf %lf %g\n",1.0/T, energy, BlockModularity3E(part,eij,ngroup,nlink,nnod), T); *\/ */
/* /\*     printf("%g %lf %lf %g\n",1.0/T, energy, BlockModularity3(part,net), T); *\/ */
/* /\*     if (energy != BlockModularity3(part,net)) *\/ */
/* /\*       printf("\n\n\n\n\n\n\nERROR!!!!\n\n\n\n\n\n\n\n\n\n"); *\/ */

/*     for ( k=0; k < (fac*nnod*ngroup); k++ ){ */

/*       ////////////////////////////////// */
/*       // Propose an individual change // */
/*       ////////////////////////////////// */

/*       target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(gsl_rng_uniform(gen) * (double)ngroup); */
/*       }while(newg == oldg); */

/*       // Calculate the change of energy */
/*       // Old energy */
/*       eold = energy; */

/*       // Move the node and actualize the eij matrix */
/*       MoveNode(nlist[target],glist[oldg],glist[newg]); */
/*       nod = nlist[target]->neig; */
/*       while(nod->next != NULL){ */
/* 	nod = nod->next; */
/* 	eij[Map2Array(oldg,(nod->ref)->inGroup,ngroup)] -= 1; */
/* 	eij[Map2Array(newg,(nod->ref)->inGroup,ngroup)] += 1; */
/*       } */

/*       // New energy */
/*       enew = BlockModularity3E(part,eij,ngroup,nlink,nnod); */

/*       // Accept the change according to Metroppolis */
/*       dE = enew - eold; */
/*       if( (dE <= 0.0) || ( gsl_rng_uniform(gen) < exp(-dE/T) ) ){ */
/* 	energy = enew; */
/*       } */
/*       else{ */
/* 	// "Actualize back" the eij matrix and move the node back */
/* 	MoveNode(nlist[target],glist[newg],glist[oldg]); */
/* 	nod = nlist[target]->neig; */
/* 	while(nod->next != NULL){ */
/* 	  nod = nod->next; */
/* 	  eij[Map2Array(oldg,(nod->ref)->inGroup,ngroup)] += 1; */
/* 	  eij[Map2Array(newg,(nod->ref)->inGroup,ngroup)] -= 1; */
/* 	} */
/*       } */

/*     } */

/*     T = T * Ts; */
/*   } */

/*   return part; */
/* } */

/* // Identifies communities using the Girvan-Newman algorithm. */
/* // Returns the partition obtained after max_lev iterations of the */
/* // split procedure */
/* struct group *GNCommunityIdentLev(struct node_gra *root, int max_lev) */
/* { */
/*   struct node_gra *root_cop = NULL; */
/*   struct node_gra *root_loc = NULL; */
/*   struct group *part; */
/*   struct group *g; */
/*   int *n1, *n2; */
/*   int rn1, rn2; */
/*   int lev = 0; */

/*   n1 = &rn1; */
/*   n2 = &rn2; */

/*   root_cop = CopyNetwork(root); */
/*   part = ClustersPartition(root_cop); */

/*   while (lev < max_lev) { */
/*     lev++; */
/*     printf("Level %d\n",lev); */

/*     g = part; */
/*     while (g->next != NULL) { */

/*       g = g->next; */
/*       root_loc = BuildNetFromGroup(g); */

/*       printf("subnetwork size = %d\n",CountNodes(root_loc)); */

/*       if (CountNodes(root_loc) > 1){ */
/* 	do{ */
/* 	  CalculateBiggestLinkBetweenness(root_loc,n1,n2); */
/* 	  RemoveLink(GetNode(*n1,root_cop),GetNode(*n2,root_cop)); */
/* 	  RemoveLink(GetNode(*n1,root_loc),GetNode(*n2,root_loc)); */
/* 	  printf("Link %d-%d removed\n",*n1+1,*n2+1); */
/* 	}while(IsClusterConnected(root_loc) == 1); */

/* 	printf("Split!\n\n"); */
/*       } */

/*       PrintPajekGraphTranslation(root); */

/*       RemoveGraph(root_loc); */
/*     } */

/*     RemovePartition(part); */
/*     part = ClustersPartition(root_cop); */
/*   } */

/*   PrintGroups(part); */
/*   MapPartToNet(part, root); */
/*   RemoveGraph(root_cop); */

/*   return part; */
/* } */

/* // Identifies communities using the Girvan-Newman algorithm. */
/* // Returns the partition obtained after max_lev iterations of the */
/* // split procedure */
/* struct group *GNCommunityIdent(struct node_gra *root) */
/* { */
/*   struct node_gra *root_cop = NULL; */
/*   struct node_gra *root_loc = NULL; */
/*   struct group *part; */
/*   struct group *best_part; */
/*   struct group *g; */
/*   int *n1, *n2; */
/*   int rn1, rn2; */
/*   int ngroup, ngroupant; */
/*   int nnod; */
/*   double max_mod; */

/*   n1 = &rn1; */
/*   n2 = &rn2; */

/*   nnod = CountNodes(root); */
/*   root_cop = CopyNetwork(root); */

/*   ngroup = ngroupant = 1; */
/*   best_part = ClustersPartition(root_cop); */
/*   max_mod = 0.0; */

/*   while ( ngroup < nnod) { */

/*     CalculateBiggestLinkBetweenness(root_cop,n1,n2); */
/*     RemoveLink(GetNode(*n1,root_cop),GetNode(*n2,root_cop)); */
/* /\*     printf("Link %d-%d removed\n",*n1+1,*n2+1); *\/ */

/*     ngroup = CountStronglyConnectedSets(root_cop); */

/*     if (ngroup != ngroupant){ // the network has been split */
/*       ngroupant = ngroup; */

/*       // Check whether the modularity has increased or not */
/*       part = ClustersPartition(root_cop); */
/*       MapPartToNet(part, root); */
/*       printf("%d %g\n",ngroup,Modularity(part)); */
/*       if (Modularity(part) > max_mod){ */

/* 	RemovePartition(best_part); */
/* 	RemovePartition(part); */
/* 	best_part = ClustersPartition(root_cop); */
/* 	MapPartToNet(best_part, root); */
/* 	max_mod = Modularity(best_part); */

/*       } */
/*     } */

/*   } */

/* /\*   PrintGroups(best_part); *\/ */
/*   RemoveGraph(root_cop); */

/*   return best_part; */
/* } */


/* PlotPajekBlockmodel(struct node_gra *net, struct group *part) */
/* { */
/*   int trans[maxim_int]; */
/*   int i; */
/*   struct group *g; */
/*   struct group *g2; */
/*   int ngroup; */
/*   FILE *outf; */
/*   int count = 0; */
/*   double scale = 0.7; */
/*   int nlink, maxnlink , minnlink; */

/*   // Initial operations */
/*   MapPartToNet(part,net); */
/*   CompressPart(part); */
/*   ngroup = CountGroups(part); */

/*   outf = fopen("blockmodel.net","w"); */
/*   fprintf(outf,"*Vertices %d\n",ngroup); */

/*   // Set the translation table and print the nodes with sizes */
/*   g = part; */
/*   while (g->next != NULL){ */
/*     g = g->next; */

/*     trans[g->label] = count++; */

/*     fprintf(outf,"%d \"%d\"      ellipse x_fact %g y_fact %g\n", */
/* 	    count, g->label+1, */
/* 	    scale * sqrt(g->size), scale * sqrt(g->size)); */
/*   } */

/*   // Calculate max and min G2G links */
/*   maxnlink = 0; */
/*   minnlink = CountNodes(net) * CountNodes(net); */

/*   g = part; */
/*   while (g->next != NULL){ */
/*     g = g->next; */
/*     g2 = g; */

/*     while (g2->next != NULL){ */
/*       g2 = g2->next; */

/*       nlink = NG2GLinks(g,g2); */

/*       if (nlink > 0){ */
/* 	if (nlink > maxnlink) maxnlink = nlink; */
/* 	if (nlink < minnlink) minnlink = nlink; */
/*       } */
/*     } */
/*   } */

/*   // Print the links with widths */
/*   fprintf(outf,"*Arcs\n"); */
/*   fprintf(outf,"*Edges\n"); */

/*   g = part; */
/*   while (g->next != NULL){ */
/*     g = g->next; */
/*     g2 = g; */

/*     while (g2->next != NULL){ */
/*       g2 = g2->next; */

/*       nlink = NG2GLinks(g,g2); */

/*       if (nlink > 0) { */
/* 	fprintf(outf,"%d %d %g w %g\n", trans[g->label]+1, */
/* 		trans[g2->label]+1, 4.9 * (double)(nlink-minnlink) / */
/* 		(double)(maxnlink-minnlink) + 0.1, */
/* 		4.9 * (double)(nlink-minnlink) / */
/* 		(double)(maxnlink-minnlink) + 0.1); */
/*       } */
/*     } */
/*   } */

/*   fclose(outf); */

/* } */


/* PlotPajekModuleRole(struct node_gra *net, struct group *gpart, struct group *rpart) */
/* { */
/*   struct group *pajekpart; */
/*   struct group *g; */
/*   struct group *g2; */
/*   struct node_lis *nod; */
/*   int gcount = 0; */
/*   int pajekcolor[max_size]; */
/*   int i; */
/*   FILE *outf; */

/*   // Create the pajek partition. In this partition, each group */
/*   // represents a module except for nodes with roles 3-4 and 6-7 in */
/*   // the catalog, which appear as single node groups. */
/*   pajekpart = CreateHeaderGroup(); */

/*   // Point the pointers in the partition to the network and */
/*   // set the inGroup field of the nodes according to the role partition */
/*   MapPartToNet(gpart,net); */
/*   MapPartToNet(rpart,net); */

/*   // Place each node in a partition */
/*   g = gpart; */
/*   while (g->next != NULL) { */
/*     g = g->next; */

/*     g2 = CreateGroup(pajekpart,gcount++); */
/*     pajekcolor[gcount-1] = g->label+1; */

/*     nod = g->nodeList; */
/*     while (nod->next != NULL) { */
/*       nod = nod->next; */

/*       if ( (nod->ref->inGroup !=2) && (nod->ref->inGroup !=3) && */
/* 	   (nod->ref->inGroup !=5) && (nod->ref->inGroup !=6)) { */
/* 	AddNodeToGroup(g2, nod->ref); */
/*       } */
/*       else { */
/* 	AddNodeToGroup(CreateGroup(pajekpart,gcount++),nod->ref); */
/* 	pajekcolor[gcount-1] = g->label+1; */
/*       } */
/*     } */
/*   } */

/*   MapPartToNet(pajekpart, net); */

/*   PlotPajekBlockmodel(net, pajekpart); */

/*   outf = fopen("blockmodel.clu","w"); */
/*   fprintf(outf,"*Vertices %d\n",gcount); */
/*   for (i=0; i<gcount; i++){ */
/*     fprintf(outf,"   %d\n",pajekcolor[i]); */
/*   } */
/*   fclose(outf); */

/*   MapPartToNet(gpart, net); */
/*   RemovePartition(pajekpart); */
/* } */

/* PlotXfigModuleRole(struct node_gra *net, struct group *gpart, struct group *rpart, double mod_pos[][2], gsl_rng *gen) */
/* { */
/*   struct group *g, *g2, *g3, *parentg; */
/*   struct group *xfigpart; */
/*   struct node_lis *nod; */
/*   struct node_gra *anode; */
/*   struct group *agroup; */
/*   struct group *agroup2; */
/*   struct group *specialg[max_size]; */
/*   int nspec; */
/*   FILE *outf; */
/*   int size; */
/*   int fac = 60; */
/*   double scale = 10000; */
/*   int xcent, ycent; */
/*   int i, j, iter, iterfac = 10000; */
/*   int Ei, Ef; */
/*   double dx, dy, dis; */
/*   int target; */
/*   int nlink; */
/*   double linefac = 0.7; */
/*   int startx, starty, endx, endy; */
/*   double theta; */
/*   int sign, active; */
/*   int role; */
/*   int linktype, color; */

/*   xfigpart = CreateHeaderGroup(); */

/*   outf = fopen("modules-roles.fig","w"); */

/*   fprintf(outf,"#FIG 3.2\nLandscape\nCenter\n"); */
/*   fprintf(outf,"Inches\nLetter\n100.00\nSingle\n-2\n1200 2\n"); */

/*   // Point the pointers in the partition to the network and */
/*   // set the inGroup field of the nodes according to the role partition */
/*   MapPartToNet(gpart,net); */
/*   MapPartToNet(rpart,net); */

/*   // Place each node in the proper xfigpart partition */
/*   g = gpart; */
/*   while (g->next != NULL) { */
/*     g = g->next; */

/*     // Create a module in the xfigpart for each module in the gpart */
/*     g2 = CreateGroup(xfigpart,g->label); */

/*     size = floor(0.5 + (double)fac * sqrt(g->size)); */
/*     xcent = floor(0.5 + scale*mod_pos[g->label][0]); */
/*     ycent = floor(0.5 + scale*mod_pos[g->label][1]); */

/*     g2->coorX = xcent; */
/*     g2->coorY = ycent; */

/*     // Plot the module as a node in xfig */
/*     fprintf(outf, */
/* 	    "1 3 0 1 0 7 50 0 20 0.000 1 0.0000 %d %d %d %d %d %d %d %d\n", */
/* 	    xcent, ycent, size, size, xcent, ycent, xcent+size, ycent); */

/*     // Look for nodes with special roles */
/*     nod = g->nodeList; */
/*     while (nod->next != NULL) { */
/*       nod = nod->next; */
/*       anode = nod->ref; */

/*       if ( anode->inGroup == 2 || anode->inGroup == 3 || */
/* 	   anode->inGroup == 5 || anode->inGroup == 6 ) { */
/* 	// create a group for this special module only */

/* 	g3 = CreateGroup(xfigpart,1000000+anode->num); */

/* 	// In the z coordinate we store the role */
/* 	g3->coorZ = anode->inGroup; */

/* 	AddNodeToGroup(g3,anode); */

/* 	// increase the size of the parent module, too */
/* 	g2->size += 1; */

/* 	// place the special node at random inside the parent module */
/* 	anode->coorX = xcent + */
/* 	  0.5 * (double)(size * (gsl_rng_uniform(gen) - 0.5)); */
/* 	anode->coorY = ycent + */
/* 	  0.5 * (double)(size * (gsl_rng_uniform(gen) - 0.5)); */

/* 	g3->coorX = anode->coorX; */
/* 	g3->coorY = anode->coorY; */
/* 	// The coordinates are now stored both in the node and */
/* 	// in the corresponding group in the xfigpart partition */
/*       } */

/*       else{ // if the node is not special, add to the main group */
/* 	AddNodeToGroup(g2,anode); */
/*       } */
/*     } */
/*   } */

/*   // Place the special nodes properly and print them */
/*   MapPartToNet(xfigpart,net); */

/*   g = xfigpart; */
/*   nspec = 0; */

/*   while (g->next != NULL) { */
/*     g = g->next; */

/*     if (g->size > 1) { */

/*       if (nspec > 0) { */

/* 	size = floor( 0.5 + fac * sqrt(parentg->size) ); */
/* 	xcent = floor( 0.5 + parentg->coorX); */
/* 	ycent = floor( 0.5 + parentg->coorY); */

/* 	iter = nspec * iterfac; */

/* 	for (i=0; i<iter; i++){ */

/* 	  target = floor(gsl_rng_uniform(gen) * (double)nspec); */
/* 	  agroup = specialg[target]; */

/* 	  do{ */
/* 	    dx = 0.01 * (double)(size * (gsl_rng_uniform(gen) - 0.5)); */
/* 	    dy = 0.01 * (double)(size * (gsl_rng_uniform(gen) - 0.5)); */
/* 	    dis = sqrt( (agroup->coorX + dx - xcent) * */
/* 			(agroup->coorX + dx - xcent) + */
/* 			(agroup->coorY + dy - ycent) * */
/* 			(agroup->coorY + dy - ycent) ); */
/* 	  }while (dis+fac > 0.9*(double)size); */

/* 	  Ei = Ef = 0; */

/* 	  // Contribution I: Hard discs with other special nodes */
/* 	  for (j=0; j<nspec; j++){ */
/* 	    agroup2 = specialg[j]; */
/* 	    // Old energy */
/* 	    dis = sqrt( (agroup->coorX - agroup2->coorX) * */
/* 			(agroup->coorX - agroup2->coorX) + */
/* 			(agroup->coorY - agroup2->coorY) * */
/* 			(agroup->coorY - agroup2->coorY) ); */
/* 	    if (dis < 2*fac) Ei += 100; */

/* 	    // New energy */
/* 	    dis = sqrt( (agroup->coorX + dx - agroup2->coorX) * */
/* 			(agroup->coorX + dx - agroup2->coorX) + */
/* 			(agroup->coorY + dy - agroup2->coorY) * */
/* 			(agroup->coorY + dy - agroup2->coorY) ); */
/* 	    if (dis < 2*fac) Ef += 100; */
/* 	  } */

/* 	  // Contribution II: Springs in links */
/* 	  g3 = xfigpart; */
/* 	  while (g3->next != NULL) { */
/* 	    g3 = g3->next; */

/* 	    if (g3 == parentg) // skip parent module */
/* 	      g3 = g3->next; */

/* 	    // Old energy */
/* 	    dis = sqrt( (agroup->coorX - g3->coorX) * */
/* 			(agroup->coorX - g3->coorX) + */
/* 			(agroup->coorY - g3->coorY) * */
/* 			(agroup->coorY - g3->coorY) ); */
/* 	    Ei += 1.e-6 * NG2GLinks(agroup,g3) * dis * dis; */

/* 	    // New energy */
/* 	    dis = sqrt( (agroup->coorX + dx - g3->coorX) * */
/* 			(agroup->coorX + dx - g3->coorX) + */
/* 			(agroup->coorY + dy - g3->coorY) * */
/* 			(agroup->coorY + dy - g3->coorY) ); */
/* 	    Ef += 1.e-6 * NG2GLinks(agroup,g3) * dis * dis; */
/* 	  } */

/* 	  // Stepest descend acceptance */
/* 	  if (Ef <= Ei){ */
/* 	    agroup->coorX += dx; */
/* 	    agroup->coorY += dy; */
/* 	  } */
/* 	} */

/* 	// Plot the special nodes */
/* 	for (j=0; j<nspec; j++){ */
/* 	  agroup2 = specialg[j]; */
/* 	  g3 = specialg[j]; */

/* 	  xcent = floor(agroup2->coorX); */
/* 	  ycent = floor(agroup2->coorY); */

/* 	  g3->coorX = (double)xcent; */
/* 	  g3->coorY = (double)ycent; */

/* 	  role = (int)floor(g3->coorZ + 0.5); */

/* 	  if ( role == 2 || role == 3 ){ */
/* 	    fprintf(outf, */
/* 		    "2 5 0 1 0 -1 38 0 -1 0.000 0 0 -1 0 0 5\n"); */
/* 	    fprintf(outf,"0 role3.gif\n"); */
/* 	    fprintf(outf,"%d %d %d %d %d %d %d %d %d %d\n", */
/* 		    xcent-89, ycent-77, xcent+89, ycent-77, */
/* 		    xcent+89, ycent+77, xcent-89, ycent+77, */
/* 		    xcent-89, ycent-77); */
/* 	  } */

/* 	  if ( role == 5 || role == 6 ){ */
/* 	    fprintf(outf, */
/* 		    "2 5 0 1 0 -1 38 0 -1 0.000 0 0 -1 0 0 5\n"); */
/* 	    fprintf(outf,"0 role6.gif\n"); */
/* 	    fprintf(outf,"%d %d %d %d %d %d %d %d %d %d\n", */
/* 		    xcent-83, ycent-83, xcent+83, ycent-83, */
/* 		    xcent+83, ycent+83, xcent-83, ycent+83, */
/* 		    xcent-83, ycent-83); */
/* 	  } */

/* 	  fprintf(outf, */
/* 		  "4 0 0 35 0 16 12 0.0000 4 135 390 %d %d %d\\001\n", */
/* 		  xcent, ycent, g3->label+1); */
/* 	  printf("Node: %d\n",g3->label+1); */
/* 	} */
/*       } // end of: if (nspec > 0) */

/*       parentg = g; // save the last g with more than one node */
/*       nspec = 0;   // reset */

/*     } // end of: if (g->size > 1) */
/*     else{ */
/*       specialg[nspec++] = g; */
/*     } */

/*   } // end of loop over the groups in the xfigpart */

/*   // THE LAST GROUP NEEDS TO BE DONE SEPARATELY!!!!!!! */
/*   if (nspec > 0) { */

/*     size = floor( 0.5 + fac * sqrt(parentg->size) ); */
/*     xcent = floor( 0.5 + parentg->coorX); */
/*     ycent = floor( 0.5 + parentg->coorY); */

/*     iter = nspec * iterfac; */

/*     for (i=0; i<iter; i++){ */

/*       target = floor(gsl_rng_uniform(gen) * (double)nspec); */
/*       agroup = specialg[target]; */

/*       do{ */
/* 	dx = 0.01 * (double)(size * (gsl_rng_uniform(gen) - 0.5)); */
/* 	dy = 0.01 * (double)(size * (gsl_rng_uniform(gen) - 0.5)); */
/* 	dis = sqrt( (agroup->coorX + dx - xcent) * */
/* 		    (agroup->coorX + dx - xcent) + */
/* 		    (agroup->coorY + dy - ycent) * */
/* 		    (agroup->coorY + dy - ycent) ); */
/*       }while (dis+fac > 0.9*(double)size); */

/*       Ei = Ef = 0; */

/*       // Contribution I: Hard discs with other special nodes */
/*       for (j=0; j<nspec; j++){ */
/* 	agroup2 = specialg[j]; */
/* 	// Old energy */
/* 	dis = sqrt( (agroup->coorX - agroup2->coorX) * */
/* 		    (agroup->coorX - agroup2->coorX) + */
/* 		    (agroup->coorY - agroup2->coorY) * */
/* 		    (agroup->coorY - agroup2->coorY) ); */
/* 	if (dis < 2*fac) Ei += 100; */

/* 	// New energy */
/* 	dis = sqrt( (agroup->coorX + dx - agroup2->coorX) * */
/* 		    (agroup->coorX + dx - agroup2->coorX) + */
/* 		    (agroup->coorY + dy - agroup2->coorY) * */
/* 		    (agroup->coorY + dy - agroup2->coorY) ); */
/* 	if (dis < 2*fac) Ef += 100; */
/*       } */

/*       // Contribution II: Springs in links */
/*       g3 = xfigpart; */
/*       while (g3->next != NULL) { */
/* 	g3 = g3->next; */

/* 	if (g3 == parentg) // skip parent module */
/* 	  g3 = g3->next; */

/* 	// Old energy */
/* 	dis = sqrt( (agroup->coorX - g3->coorX) * */
/* 		    (agroup->coorX - g3->coorX) + */
/* 		    (agroup->coorY - g3->coorY) * */
/* 		    (agroup->coorY - g3->coorY) ); */
/* 	Ei += 1.e-6 * NG2GLinks(agroup,g3) * dis * dis; */

/* 	// New energy */
/* 	dis = sqrt( (agroup->coorX + dx - g3->coorX) * */
/* 		    (agroup->coorX + dx - g3->coorX) + */
/* 		    (agroup->coorY + dy - g3->coorY) * */
/* 		    (agroup->coorY + dy - g3->coorY) ); */
/* 	Ef += 1.e-6 * NG2GLinks(agroup,g3) * dis * dis; */
/*       } */

/*       // Stepest descend acceptance */
/*       if (Ef <= Ei){ */
/* 	agroup->coorX += dx; */
/* 	agroup->coorY += dy; */
/*       } */
/*     } */

/*     // Plot the special nodes */
/*     for (j=0; j<nspec; j++){ */
/*       agroup2 = specialg[j]; */
/*       g3 = specialg[j]; */

/*       xcent = floor(agroup2->coorX); */
/*       ycent = floor(agroup2->coorY); */

/*       g3->coorX = (double)xcent; */
/*       g3->coorY = (double)ycent; */

/*       role = (int)floor(g3->coorZ + 0.5); */

/*       if ( role == 2 || role == 3 ){ */
/* 	fprintf(outf, */
/* 		"2 5 0 1 0 -1 38 0 -1 0.000 0 0 -1 0 0 5\n"); */
/* 	fprintf(outf,"0 role3.gif\n"); */
/* 	fprintf(outf,"%d %d %d %d %d %d %d %d %d %d\n", */
/* 		xcent-89, ycent-77, xcent+89, ycent-77, */
/* 		xcent+89, ycent+77, xcent-89, ycent+77, */
/* 		xcent-89, ycent-77); */
/*       } */

/*       if ( role == 5 || role == 6 ){ */
/* 	fprintf(outf, */
/* 		"2 5 0 1 0 -1 38 0 -1 0.000 0 0 -1 0 0 5\n"); */
/* 	fprintf(outf,"0 role6.gif\n"); */
/* 	fprintf(outf,"%d %d %d %d %d %d %d %d %d %d\n", */
/* 		xcent-83, ycent-83, xcent+83, ycent-83, */
/* 		xcent+83, ycent+83, xcent-83, ycent+83, */
/* 		xcent-83, ycent-83); */
/*       } */

/*       fprintf(outf, */
/* 	      "4 0 0 35 0 16 12 0.0000 4 135 390 %d %d %d\\001\n", */
/* 	      xcent, ycent, g3->label+1); */
/*       printf("Node: %d\n",g3->label+1); */
/*     } */
/*   } // end of: if (nspec > 0) */



/*   // Plot the links using the xfigpart partition */

/*   printf("--------------\n"); */
/*   PrintGroups(xfigpart); */
/*   printf("--------------\n"); */

/*   g2 = xfigpart; */

/*   while (g2->next != NULL){ */
/*     g2 = g2->next; */

/*     active = 0; // to avoid printing the link from a module */
/*                 // to its special nodes */

/*     g3 = g2; */
/*     while (g3->next != NULL){ */
/*       g3 = g3->next; */

/*       if (g3->size > 1 ) active = 1; */

/*       // determine the origin of the line */
/*       if (active == 1){ */

/* 	if (g2->size == 1) { */
/* 	  startx = g2->coorX; */
/* 	  starty = g2->coorY; */
/* 	} */
/* 	else{ */
/* 	  theta = atan( (g3->coorY - g2->coorY) / */
/* 			(g3->coorX - g2->coorX) ); */
/* 	  if (g3->coorX > g2->coorX) */
/* 	    sign = +1; */
/* 	  else */
/* 	    sign = -1; */

/* 	  startx = floor(0.5 + g2->coorX + sign * */
/* 			 cos(theta) * fac * sqrt(g2->size)); */
/* 	  starty = floor(0.5 + g2->coorY + sign * */
/* 			 sin(theta) * fac * sqrt(g2->size)); */
/* 	} */

/*       // determine the destination of the line */
/* 	if (g3->size == 1) { */
/* 	  endx = g3->coorX; */
/* 	  endy = g3->coorY; */
/* 	} */
/* 	else{ */
/* 	  theta = atan( (g3->coorY - g2->coorY) / */
/* 			(g3->coorX - g2->coorX) ); */
/* 	  if (g3->coorX > g2->coorX) */
/* 	    sign = -1; */
/* 	  else */
/* 	    sign = +1; */

/* 	  endx = floor(0.5 + g3->coorX + sign * */
/* 		       cos(theta) * fac * sqrt(g3->size)); */
/* 	  endy = floor(0.5 + g3->coorY + sign * */
/* 		       sin(theta) * fac * sqrt(g3->size)); */
/* 	} */

/* 	// count the number of links */
/* 	nlink = NG2GLinks(g2,g3); */

/* 	linktype = 0; */
/* 	if (g2->size == 1){ */
/* 	  linktype++; */
/* 	} */
/* 	if (g3->size == 1){ */
/* 	  linktype++; */
/* 	} */
/* 	if (linktype == 0) */
/* 	  color = 4; */
/* 	else if (linktype == 1) */
/* 	  color = 1; */
/* 	else */
/* 	  color = 0; */

/* 	// print the link */
/* 	if ( nlink > 0 ) { */
/* 	  printf("%d (%d) - %d (%d): %d links\n", */
/* 		 g2->label+1,g2->size,g3->label+1,g3->size,nlink); */

/* 	  fprintf(outf,"2 1 0 %d %d 7 45 0 -1 0.000 0 0 -1 0 0 2\n", */
/* 		  (int)(ceil)(linefac * (double)nlink), color); */
/* 	  fprintf(outf,"\t %d %d %d %d\n", startx, starty, endx, endy); */
/* 	} */
/*       } //end of: if (active == 1) */
/*     } // end of g3 loop */
/*   } // end of g2 loop */

/*   fclose(outf); */

/*   RemovePartition(xfigpart); */
/* } */

/* // Randomizes the links in a network keeping the role of each */
/* // node. This is achieved by randomizing links between modules i and j */
/* // with each other---but not with other links---for all i and j */
/* // (including j=i). */
/* struct node_gra *RandomizeNetworkRolePreserve(struct node_gra *net, struct group *part, double times, gsl_rng *gen) */
/* { */
/*   int i, j; */
/*   int nlink; */
/*   int target1, target2; */
/*   struct node_gra *p = NULL; */
/*   struct node_lis *l = NULL; */
/*   struct node_gra *n1, *n2, *n3, *n4; */
/*   struct node_gra **ori, **des; */
/*   int coun = 0; */
/*   int niter; */
/*   struct group *g1, *g2; */
/*   int randomizable; */

/*   // ------------------------------ */
/*   // Loop over all pairs of modules */
/*   // ------------------------------ */
/*   g1 = part; */
/*   while (g1->next != NULL){ */
/*     g2 = g1; */
/*     g1 = g1->next; */
/*     while (g2->next != NULL){ */
/*       g2 = g2->next; */

/* /\*       printf("\tModules %d and %d\n", g1->label+1, g2->label+1); *\/ */

/*       // Build the link lists (one for link origins and one for ends) */
/*       nlink = NG2GLinks(g1, g2); */
/*       if (g1 == g2) */
/* 	nlink /= 2; */
/*       niter =  ceil(times * (double)nlink + 0.5 + gsl_rng_uniform(gen)); */
/*       // The +0.5+ran2 makes the number of iterations sometimes even */
/*       // and sometimes odd, for cases in which thre may be only one */
/*       // pair of randomizable links. */

/*       ori = (struct node_gra **) */
/* 	calloc(nlink,sizeof(struct node_gra *)); */
/*       des = (struct node_gra **) */
/* 	calloc(nlink,sizeof(struct node_gra *)); */

/*       // Find the links and add them to the lists ori and des */
/*       p = net; */
/*       coun = 0; */
/*       while(p->next !=  NULL){ */
/* 	p = p->next; */

/* 	l = p->neig; */
/* 	while(l->next !=  NULL){ */
/* 	  l = l->next; */

/* 	  if ( ( (g1 == g2) && (p->num > l->node) && */
/* 		 (p->inGroup == g1->label) && */
/* 		 (l->ref->inGroup == g2->label) */
/* 		 ) || */
/* 	       ( (g1 != g2) &&  */
/* 		 (p->inGroup == g1->label) && */
/* 		 (l->ref->inGroup == g2->label) ) */
/* 	       ) { */
/* 	    ori[coun] = p; */
/* 	    des[coun] = l->ref; */
/* 	    coun++; */
/* 	  } */
/* 	} */
/*       } */

/* /\*       printf("%d %d %d %d\n", g1->label+1, g2->label+1, nlink, coun); *\/ */

/*       if(coun !=  nlink) */
/* 	printf("Error in RandomizeNetworkPreserveRole: coun !=  nlink!!\n"); */

/*       // Check if there is at least a pair of randomizable links */
/*       randomizable = 0; */
/*       for ( i=0; i<nlink; i++ ) { */
/* 	for ( j=i+1; j<nlink; j++ ) { */
/* 	  n1 = ori[i]; */
/* 	  n2 = des[i]; */
/* 	  n3 = ori[j]; */
/* 	  n4 = des[j]; */
/* 	  if ( n3 != n1 && n3 != n2 && n4 != n1 && n4 != n2 && */
/* 	       IsThereLink(n1,n4->num) != 1 && */
/* 	       IsThereLink(n2,n3->num) != 1 ) { */
/* 	    randomizable = 1; */
/* 	    i = j = nlink; */
/* 	  } */
/* 	} */
/*       } */

/*       if ( randomizable == 1 ) { */
/* 	// Randomize the links */
/* 	for(i = 0; i<niter; i++){ */

/* 	  // select the 4 different nodes */
/* 	  do{ */
/* 	    target1 = floor(gsl_rng_uniform(gen) * (double)nlink); */
/* 	    n1 = ori[target1]; */
/* 	    n2 = des[target1]; */

/* /\* 	    printf("%d-%d\n", n1->num+1, n2->num+1); *\/ */

/* 	    target2 = floor(gsl_rng_uniform(gen) * (double)nlink); */
/* 	    if ( (g1 == g2) && (gsl_rng_uniform(gen) < 0.5) ) { */
/* 	      n4 = ori[target2]; */
/* 	      n3 = des[target2]; */
/* 	    } */
/* 	    else { */
/* 	      n3 = ori[target2]; */
/* 	      n4 = des[target2]; */
/* 	    } */
/* 	  }while(n3 == n1 || n3 == n2 || n4 == n1 || n4 == n2 || */
/* 		 IsThereLink(n1,n4->num) == 1 || */
/* 		 IsThereLink(n2,n3->num) == 1); */

/* /\* 	  printf("%d-%d and %d-%d\n", *\/ */
/* /\* 		 n1->num+1, n2->num+1, n3->num+1, n4->num+1); *\/ */

/* 	  // switch the link */
/* 	  RemoveLink(n1,n2); */
/* 	  RemoveLink(n3,n4); */
/* 	  AddAdjacencyFull(n1,n4,1); */
/* 	  AddAdjacencyFull(n4,n1,1); */
/* 	  AddAdjacencyFull(n3,n2,1); */
/* 	  AddAdjacencyFull(n2,n3,1); */

/* 	  ori[target1] = n1; */
/* 	  des[target1] = n4; */
/* 	  ori[target2] = n3; */
/* 	  des[target2] = n2; */
/* /\* 	  printf("%d-%d %d-%d --> %d-%d %d-%d\n", *\/ */
/* /\* 		 n1->num+1, n2->num+1, n3->num+1, n4->num+1, *\/ */
/* /\* 		 n1->num+1, n4->num+1, n3->num+1, n2->num+1); *\/ */
/* 	} */
/*       } */
/* /\*       else *\/ */
/* /\* 	printf("\t\tnothing to randomize\n"); *\/ */

/*       free(ori); */
/*       free(des); */

/*     }  // End the loop over pairs of modules */
/*   } */

/*   return net; */
/* } */


/* // Randomizes the links in a network keeping the role of each */
/* // node. This is achieved by randomizing links between modules i and j */
/* // with each other---but not with other links---for all i and j */
/* // (including j=i). Contrary to RandomizeNetworkRolePreserve, it does */
/* // not allow switches that break the network (assumed to be connected) */
/* // into multiple disconnected components. */
/* struct node_gra *RandomizeNetworkOneComponentRolePreserve(struct node_gra *net, */
/* 							  struct group *part, */
/* 							  double times, */
/* 							  gsl_rng *gen) */
/* { */
/*   int i, j; */
/*   int nlink; */
/*   int target1, target2; */
/*   struct node_gra *p = NULL; */
/*   struct node_lis *l = NULL; */
/*   struct node_gra *n1, *n2, *n3, *n4; */
/*   struct node_gra **ori, **des; */
/*   int coun = 0; */
/*   int niter; */
/*   struct group *g1, *g2; */
/*   int randomizable; */
/*   int next; */

/*   // ------------------------------ */
/*   // Loop over all pairs of modules */
/*   // ------------------------------ */
/*   g1 = part; */
/*   while (g1->next != NULL){ */
/*     g2 = g1; */
/*     g1 = g1->next; */
/*     while (g2->next != NULL){ */
/*       g2 = g2->next; */

/* /\*       printf("\tModules %d and %d\n", g1->label+1, g2->label+1); *\/ */

/*       // Build the link lists (one for link origins and one for ends) */
/*       nlink = NG2GLinks(g1, g2); */
/*       if (g1 == g2) */
/* 	nlink /= 2; */
/*       niter =  ceil(times * (double)nlink + 0.5 + gsl_rng_uniform(gen)); */
/*       // The +0.5+ran2 makes the number of iterations sometimes even */
/*       // and sometimes odd, for cases in which thre may be only one */
/*       // pair of randomizable links. */

/*       ori = (struct node_gra **) */
/* 	calloc(nlink,sizeof(struct node_gra *)); */
/*       des = (struct node_gra **) */
/* 	calloc(nlink,sizeof(struct node_gra *)); */

/*       // Find the links and add them to the lists ori and des */
/*       p = net; */
/*       coun = 0; */
/*       while(p->next !=  NULL){ */
/* 	p = p->next; */

/* 	l = p->neig; */
/* 	while(l->next !=  NULL){ */
/* 	  l = l->next; */

/* 	  if ( ( (g1 == g2) && (p->num > l->node) && */
/* 		 (p->inGroup == g1->label) && */
/* 		 (l->ref->inGroup == g2->label) */
/* 		 ) || */
/* 	       ( (g1 != g2) &&  */
/* 		 (p->inGroup == g1->label) && */
/* 		 (l->ref->inGroup == g2->label) ) */
/* 	       ) { */
/* 	    ori[coun] = p; */
/* 	    des[coun] = l->ref; */
/* 	    coun++; */
/* 	  } */
/* 	} */
/*       } */

/* /\*       printf("%d %d %d %d\n", g1->label+1, g2->label+1, nlink, coun); *\/ */

/*       if(coun !=  nlink) */
/* 	printf("Error in RandomizeNetworkPreserveRole: coun !=  nlink!!\n"); */

/*       // Check if there is at least a pair of randomizable links */
/*       randomizable = 0; */
/*       for ( i=0; i<nlink; i++ ) { */
/* 	for ( j=i+1; j<nlink; j++ ) { */
/* 	  n1 = ori[i]; */
/* 	  n2 = des[i]; */
/* 	  n3 = ori[j]; */
/* 	  n4 = des[j]; */
/* 	  if ( n3 != n1 && n3 != n2 && n4 != n1 && n4 != n2 && */
/* 	       IsThereLink(n1,n4->num) != 1 && */
/* 	       IsThereLink(n2,n3->num) != 1 ) { */

/* 	    // Check if this switch will leave the network connected */
/* 	    RemoveLink(n1, n2); */
/* 	    RemoveLink(n3, n4); */
/* 	    AddAdjacencyFull(n1, n4, 1); */
/* 	    AddAdjacencyFull(n4, n1, 1); */
/* 	    AddAdjacencyFull(n3, n2, 1); */
/* 	    AddAdjacencyFull(n2, n3, 1); */

/* 	    if (IsGraphConnected(net) == 1) { // Is the network still */
/* 					      // connected? */
/* 	      randomizable = 1;  */
/* 	      i = j = nlink;     // done */
/* 	      fprintf(stderr, "%d-%d %d-%d\n", */
/* 		      n1->num+1, n2->num+1, n3->num+1, n4->num+1); */
/* 	    } */

/* 	    // Undo the switch */
/* 	    RemoveLink(n1, n4); */
/* 	    RemoveLink(n2, n3); */
/* 	    AddAdjacencyFull(n1, n2, 1); */
/* 	    AddAdjacencyFull(n2, n1, 1); */
/* 	    AddAdjacencyFull(n3, n4, 1); */
/* 	    AddAdjacencyFull(n4, n3, 1); */
/* 	  } */
/* 	} */
/*       } */

/*       if ( randomizable == 1 ) { */
/* 	// Randomize the links */
/* 	for(i = 0; i<niter; i++){ */

/* 	  do { // loop to make sure that the switch does not break the */
/* 	       // network into disconnected components */
/* 	    next = 0; */

/* 	    // select the 4 different nodes */
/* 	    do{ */
/* 	      target1 = floor(gsl_rng_uniform(gen) * (double)nlink); */
/* 	      n1 = ori[target1]; */
/* 	      n2 = des[target1]; */

/* /\* 	    printf("%d-%d\n", n1->num+1, n2->num+1); *\/ */

/* 	      target2 = floor(gsl_rng_uniform(gen) * (double)nlink); */
/* 	      if ( (g1 == g2) && (gsl_rng_uniform(gen) < 0.5) ) { */
/* 		n4 = ori[target2]; */
/* 		n3 = des[target2]; */
/* 	      } */
/* 	      else { */
/* 		n3 = ori[target2]; */
/* 		n4 = des[target2]; */
/* 	      } */
/* 	    }while(n3 == n1 || n3 == n2 || n4 == n1 || n4 == n2 || */
/* 		   IsThereLink(n1,n4->num) == 1 || */
/* 		   IsThereLink(n2,n3->num) == 1); */

/* /\* 	  printf("%d-%d and %d-%d\n", *\/ */
/* /\* 		 n1->num+1, n2->num+1, n3->num+1, n4->num+1); *\/ */

/* 	    // switch the link */
/* 	    RemoveLink(n1,n2); */
/* 	    RemoveLink(n3,n4); */
/* 	    AddAdjacencyFull(n1,n4,1); */
/* 	    AddAdjacencyFull(n4,n1,1); */
/* 	    AddAdjacencyFull(n3,n2,1); */
/* 	    AddAdjacencyFull(n2,n3,1); */

/* 	    if (IsGraphConnected(net) == 1) { // Is the network still */
/* 					      // connected? */
/* 	      ori[target1] = n1; */
/* 	      des[target1] = n4; */
/* 	      ori[target2] = n3; */
/* 	      des[target2] = n2; */
/* 	      next = 1; // We have completed the switch */
/* 	    } */
/* 	    else { // Network is not connected any more: undo the link */
/* 	      RemoveLink(n1, n4); */
/* 	      RemoveLink(n2, n3); */
/* 	      AddAdjacencyFull(n1, n2, 1); */
/* 	      AddAdjacencyFull(n2, n1, 1); */
/* 	      AddAdjacencyFull(n3, n4, 1); */
/* 	      AddAdjacencyFull(n4, n3, 1); */
/* 	    } */
/* 	  } while (next == 0); // do not let go until the link has */
/* 			       // been switched */

/* /\* 	  printf("%d-%d %d-%d --> %d-%d %d-%d\n", *\/ */
/* /\* 		 n1->num+1, n2->num+1, n3->num+1, n4->num+1, *\/ */
/* /\* 		 n1->num+1, n4->num+1, n3->num+1, n2->num+1); *\/ */
/* 	} */
/*       } */
/* /\*       else *\/ */
/* /\* 	printf("\t\tnothing to randomize\n"); *\/ */

/*       free(ori); */
/*       free(des); */

/*     }  // End the loop over pairs of modules */
/*   } */

/*   return net; */
/* } */




/* /\* */
/*   Builds a network in which each node is a module of the given */
/*   partition. Links in the blockmodel network are weighted by the */
/*   number of links between the corresponding modules. */
/* *\/ */
/* struct node_gra *BuildBlockmodelNetwork(struct group *part) */
/* { */
/*   struct group *g1=NULL, *g2=NULL; */
/*   struct node_gra *net=NULL, *nlist[maxim_int], *last=NULL; */
/*   int weight; */

/*   // Create the nodes (one for each non-empty module) */
/*   last = net = CreateHeaderGraph(); */
/*   g1 = part; */
/*   while ((g1 = g1->next) != NULL) { */
/*     if (g1->size > 0) { */
/*       last = nlist[g1->label] = CreateNodeGraph(last, g1->label, 0,0); */
/*       printf("Adding %d\n", g1->label+1); */
/*     } */
/*   } */
/*   printf("done\n"); */

/*   // Create the links */
/*   g1 = part; */
/*   while ((g1 = g1->next) != NULL) { */
/*     g2 = g1; */
/*     while ((g2 = g2->next) != NULL) { */
/*       weight = NG2GLinks(g1, g2); */
/*       if (weight > 0) { */
/* 	AddWeightedAdjacencyFull(nlist[g1->label], */
/* 				 nlist[g2->label], */
/* 				 1, (double)weight); */
/* 	AddWeightedAdjacencyFull(nlist[g2->label], */
/* 				 nlist[g1->label], */
/* 				 1, (double)weight); */
/*       } */
/*     } */
/*   } */

/*   return net; */
/* } */


/* /\* */
/*   Returns a randomly selected node from the specified module */
/* *\/ */
/* struct node_gra *SelectRandomNodeInGroup(struct group *mod, gsl_rng *gen) */
/* { */
/*   int i; */
/*   struct node_lis *node = mod->nodeList->next; */
/*   int target = floor(gsl_rng_uniform(gen) * mod->size); */

/*   for (i=0; i<target; i++) { */
/*     node = node->next; */
/*   } */

/*   return node->ref; */
/* } */


/* /\* */
/*   Returns 1 if the node is in the group */
/* *\/ */
/* int IsNodeInGroup(struct node_gra *node, struct group *mod) */
/* { */
/*   struct node_lis *p = mod->nodeList; */

/*   while ((p = p->next) != NULL) */
/*     if (p->ref == node) */
/*       return 1; */

/*   return 0; */
/* } */













/* // Creates a Pajek file with the nodes of a module (and its neighbors) */
/* // only */
/* PrintPajekOneModule(struct node_gra *net, */
/* 		    struct group *mod, struct group *part, */
/* 		    struct group *part2,  */
/* 		    char netF[], char parF[], char par2F[]) */
/* { */
/*   struct node_gra *p = net; */
/*   struct node_lis *li; */
/*   FILE *fit1, *fit2, *fit3; */
/*   int trans[maxim_int]; */
/*   int i, co; */
/*   int nnod = 0; */

/*   for(i=0; i<maxim_int; i++){ */
/*     trans[i] = 0; */
/*   } */

/*   // Count the nodes */
/*   p = net; */
/*   while (p->next != NULL) { */
/*     p = p->next; */
/*     if ( (p->inGroup == mod->label) || */
/* 	 (NLinksToGroup(p, mod) > 0) ) */
/*       nnod ++; */
/*   } */

/*   // Open the files  */
/*   fit1 = fopen(netF, "w"); */
/*   fit2 = fopen(parF, "w"); */
/*   if (part2 != NULL) */
/*     fit3 = fopen(par2F, "w"); */

/*   fprintf(fit1,"*Vertices %d\n", nnod); */
/*   fprintf(fit2,"*Vertices %d\n", nnod); */
/*   if (part2 != NULL) */
/*     fprintf(fit3,"*Vertices %d\n", nnod); */

/*   // Write the files */
/*   co = 1; */
/*   p = net; */
/*   while (p->next != NULL) { */
/*     p = p->next; */

/*     // Print the node if it is in the module or directly connected to */
/*     // the module */
/*     if ( (p->inGroup == mod->label) || */
/* 	 (NLinksToGroup(p, mod) > 0) ) { */
/*       trans[p->num] = co++; */
/*       // Print the node */
/*       fprintf(fit1,"%d \"%d\"           %lf %lf %lf\n", */
/* 	      trans[p->num], p->num+1, p->coorX, p->coorY, p->coorZ); */
/*       // Print if the node is in the module (1) or directly connected */
/*       // to the module (2) */
/*       if (p->inGroup == mod->label) */
/* 	fprintf(fit2,"   1\n"); */
/*       else */
/* 	fprintf(fit2,"   2\n"); */
/*       // Print the additional partition, if any */
/*       if (part2 != NULL) { */
/* 	MapPartToNet(part2, net); */
/* 	fprintf(fit3,"   %d\n", p->inGroup+1); */
/* 	MapPartToNet(part, net); */
/*       } */
/*     } */
/*   } */

/*   fprintf(fit1, "*Arcs\n"); */
/*   fprintf(fit1, "*Edges\n"); */
/*   p = net; */
/*   while (p->next != NULL) { */
/*     p = p->next; */
/*     if (trans[p->num] != 0) { */
/*       li = p->neig; */
/*       while (li->next != NULL) { */
/* 	li = li->next; */
/* 	if(trans[li->node] != 0) */
/* 	  fprintf(fit1,"%d   %d   1\n",trans[p->num],trans[li->node]); */
/*       } */
/*     } */
/*   } */

/*   fclose(fit1); */
/*   fclose(fit2); */
/*   if (part2 != NULL) */
/*     fclose(fit3); */
/* } */























/* struct group *ThermalNetworkSplitWeight(struct group *targ,double Ti, double Tf, gsl_rng *gen) */
/* { */
/*   struct group *glist[2]; */
/*   struct group *split = NULL; */
/*   struct node_gra *nlist[max_size]; */
/*   struct node_gra *net = NULL; */
/*   struct node_gra *p = NULL; */
/*   struct node_p; */
/*   int nnod = 0; */
/*   int i; */
/*   int des; */
/*   int target,oldg,newg; */
/*   double innew,inold,nlink; */
/*   double totallinks=0.0; */
/*   double dE=0.0, energy=0.0; */
/*   double T, Ts = 0.95; */

/*   for(i=0; i<max_size; i++){ */
/*     nlist[i] = NULL; */
/*   } */
/*   glist[0] = NULL; */
/*   glist[1] = NULL; */

/*   printf("Don't you want to add the Perc thing?!\n"); */

/*   // Build a network from the nodes in the target group */
/*   net = BuildNetFromGroup(targ); */

/*   // Create the groups */
/*   split = CreateHeaderGroup(); */
/*   glist[0] = CreateGroup(split,0); */
/*   glist[1] = CreateGroup(split,1); */

/*   // Randomly assign the nodes to the groups */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     nlist[nnod] = p; */
/*     totallinks += NodeDegreeWeight(p); */
/*     nnod++; */

/*     des = floor(gsl_rng_uniform(gen)*2.0); */
/*     AddNodeToGroup(glist[des],p); */
/*   } */

/*   totallinks /= 2.0; */

/*   // Do the SA to "optimize" the splitting */
/*   T = Ti; */
/*   while( T > Tf){ */

/*     for (i=0; i< nnod; i++){ */
/*       target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       if(oldg == 0) */
/* 	newg = 1; */
/*       else */
/* 	newg = 0; */

/*       // Calculate the change of energy */
/*       inold = NodeDegreeWeightInGroup(nlist[target],glist[oldg]); */
/*       innew = NodeDegreeWeightInGroup(nlist[target],glist[newg]); */
/*       nlink = NodeDegreeWeight(nlist[target]); */

/*       dE = 0.0; */

/*       dE -= (double)(2 * glist[oldg]->inlinksW) / */
/* 	(double)totallinks - */
/* 	(double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) * */
/* 	(double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE -= (double)(2 * glist[newg]->inlinksW) / */
/* 	(double)totallinks - */
/* 	(double)(glist[newg]->totlinksW+glist[newg]->inlinksW) * */
/* 	(double)(glist[newg]->totlinksW+glist[newg]->inlinksW) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE += (double)(2*glist[oldg]->inlinksW - 2*inold) / */
/* 	(double)totallinks - */
/* 	(double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW - */
/* 		 nlink ) * */
/* 	(double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW - */
/* 		 nlink ) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE += (double)(2*glist[newg]->inlinksW + 2*innew) / */
/* 	(double)totallinks - */
/* 	(double)(glist[newg]->totlinksW + glist[newg]->inlinksW + */
/* 		 nlink ) * */
/* 	(double)(glist[newg]->totlinksW + glist[newg]->inlinksW + */
/* 		 nlink ) / */
/* 	((double)totallinks * (double)totallinks); */

/*       // Accept the change according to the Boltzman factor */
/*       if( (dE >= 0.0) || (gsl_rng_uniform(gen) < exp(dE/T)) ){ */
/* 	MoveNode(nlist[target],glist[oldg],glist[newg]); */
/* 	energy += dE; */
/*       } */
/*     } */

/*     T = T * Ts; */
/*   } */

/*   RemoveGraph(net); */

/*   return split; */
/* } */


/* struct group *ThermalPercNetworkSplitWeight(struct group *targ, */
/* 					    double Ti, double Tf, */
/* 					    gsl_rng *gen) */
/* { */
/*   struct group *glist[2]; */
/*   struct group *split = NULL; */
/*   struct node_gra *nlist[max_size]; */
/*   struct node_gra *net = NULL; */
/*   struct node_gra *p = NULL; */
/*   struct node_p; */
/*   int nnod = 0; */
/*   int i; */
/*   int des; */
/*   int target,oldg,newg; */
/*   double innew,inold,nlink; */
/*   double totallinks=0.0; */
/*   double dE=0.0, energy=0.0; */
/*   double T, Ts = 0.95; */
/*   int ngroups, g1, g2; */
/*   double prob = 0.5; */

/*   for(i=0; i<max_size; i++){ */
/*     nlist[i] = NULL; */
/*   } */
/*   glist[0] = NULL; */
/*   glist[1] = NULL; */

/*   // Build a network from the nodes in the target group */
/*   net = BuildNetFromGroup(targ); */

/*   // Check if the network is connected */
/*   split = ClustersPartition(net); */
/*   ngroups = CountGroups(split); */

/*   if ( ngroups > 1 && gsl_rng_uniform(gen) < prob) { // Network is not */
/* 						   // connected */

/*     // Merge groups randomly until only two are left */
/*     while (ngroups > 2) { */
/*       // Select two random groups */
/*       g1 = ceil(gsl_rng_uniform(gen)* (double)ngroups); */
/*       do { */
/* 	g2 = ceil(gsl_rng_uniform(gen)* (double)ngroups); */
/*       } while (g2 == g1); */

/*       glist[0] = split; */
/*       for(i=0; i<g1; i++) */
/* 	glist[0] = glist[0]->next; */
/*       glist[1] = split; */
/*       for(i=0; i<g2; i++) */
/* 	glist[1] = glist[1]->next; */

/*       // Merge */
/*       MergeGroups(glist[0], glist[1]); */
/*       split = CompressPart(split); */
/*       ngroups--; */
/*     } */
/*   } */

/*   else { // Network IS connected */
/*     // Remove SCS partition */
/*     RemovePartition(split); */
/*     ResetNetGroup(net); */

/*     // Create the groups */
/*     split = CreateHeaderGroup(); */
/*     glist[0] = CreateGroup(split,0); */
/*     glist[1] = CreateGroup(split,1); */

/*     // Randomly assign the nodes to the groups */
/*     p = net; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       nlist[nnod] = p; */
/*       totallinks += NodeDegreeWeight(p); */
/*       nnod++; */

/*       des = floor(gsl_rng_uniform(gen)*2.0); */
/*       AddNodeToGroup(glist[des],p); */
/*     } */

/*     totallinks /= 2.0; */

/*     // Do the SA to "optimize" the splitting */
/*     if ( totallinks > 0 ) { */
/*       T = Ti; */
/*       while( T > Tf){ */

/* 	for (i=0; i< nnod; i++){ */
/* 	  target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/* 	  oldg = nlist[target]->inGroup; */
/* 	  if(oldg == 0) */
/* 	    newg = 1; */
/* 	  else */
/* 	    newg = 0; */

/* 	  // Calculate the change of energy */
/* 	  inold = NodeDegreeWeightInGroup(nlist[target],glist[oldg]); */
/* 	  innew = NodeDegreeWeightInGroup(nlist[target],glist[newg]); */
/* 	  nlink = NodeDegreeWeight(nlist[target]); */

/* 	  dE = 0.0; */

/* 	  dE -= (double)(2 * glist[oldg]->inlinksW) / */
/* 	    (double)totallinks - */
/* 	    (double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) * */
/* 	    (double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) / */
/* 	    ((double)totallinks * (double)totallinks); */

/* 	  dE -= (double)(2 * glist[newg]->inlinksW) / */
/* 	    (double)totallinks - */
/* 	    (double)(glist[newg]->totlinksW+glist[newg]->inlinksW) * */
/* 	    (double)(glist[newg]->totlinksW+glist[newg]->inlinksW) / */
/* 	    ((double)totallinks * (double)totallinks); */

/* 	  dE += (double)(2*glist[oldg]->inlinksW - 2*inold) / */
/* 	    (double)totallinks - */
/* 	    (double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW - */
/* 		     nlink ) * */
/* 	    (double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW - */
/* 		     nlink ) / */
/* 	    ((double)totallinks * (double)totallinks); */

/* 	  dE += (double)(2*glist[newg]->inlinksW + 2*innew) / */
/* 	    (double)totallinks - */
/* 	    (double)(glist[newg]->totlinksW + glist[newg]->inlinksW + */
/* 		     nlink ) * */
/* 	    (double)(glist[newg]->totlinksW + glist[newg]->inlinksW + */
/* 		     nlink ) / */
/* 	    ((double)totallinks * (double)totallinks); */

/* 	  // Accept the change according to the Boltzman factor */
/* 	  if( (dE >= 0.0) || (gsl_rng_uniform(gen) < exp(dE/T)) ){ */
/* 	    MoveNode(nlist[target],glist[oldg],glist[newg]); */
/* 	    energy += dE; */
/* 	  } */
/* 	} */

/* 	T = T * Ts; */
/*       } // End of temperature loop */
/*     } // End if totallinks > 0 */
/*   } */

/*   RemoveGraph(net); */
/*   return split; */
/* } */














/* // merge = 0 => No group merging */
/* struct group *SACommunityIdentWeight(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int merge, gsl_rng *gen) */
/* { */
/*   int i; */
/*   struct group *part = NULL; */
/*   struct group *split = NULL, *g = NULL; */
/*   struct group *glist[max_size]; */
/*   struct node_gra *nlist[max_size]; */
/*   struct node_gra *p; */
/*   struct node_lis *nod; */
/*   int target,empty; */
/*   int newg,oldg; */
/*   int nnod; */
/*   double totallinks = 0.0; */
/*   double innew,inold,nlink; */
/*   double energy = 0.0, dE; */
/*   double T; */
/*   int g1,g2; */
/*   double energyant = 0.0; */
/*   int count = 0, limit = 25; // to stop the search if the energy */
/*                           // does not change */
/*   int cicle1,cicle2; */
/*   int trans[maxim_int]; */


/*   // Create the groups and assign each node to one group */
/*   nnod = CountNodes(net); */
/*   part = CreateHeaderGroup(); */
/*   p = net->next; */
/*   ResetNetGroup(net); // All nodes reset to group -1 */

/*   nlist[0] = p; */
/*   trans[p->num] = 0; */
/*   glist[0] = CreateGroup(part,0); */
/*   AddNodeToGroup(glist[0],p); */
/*   totallinks += NodeDegreeWeight(p); */

/*   for( i=1; i<nnod; i++ ) { */
/*     p = p->next; */

/*     nlist[i] = p; */
/*     trans[p->num] = i; */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*     AddNodeToGroup(glist[i],p); */
/*     totallinks += NodeDegreeWeight(p); */
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

/*   // Do the simulated annealing */
/*   T = Ti; */
/*   energy = ModularityWeight(part); */

/*   while( T > Tf && count < limit){ */

/*     if (fabs(energy - energyant) / fabs(energy) < 1.0e-12) */
/*       count++; */
/*     else{ */
/*       energyant = energy; */
/*       count = 0; */
/*     } */

/* /\*     PrintGroups(part); *\/ */
/* /\*     printf("%g %lf %lf %g %d\n",1.0/T, energy, ModularityWeight(part), T, CountNonEmptyGroups(part)); *\/ */
/*     printf("%g %lf %g %d\n",1.0/T, energy, T, CountNonEmptyGroups(part)); */

/*     if (merge == 1) { */

/*       for ( i=0; i < cicle2; i++ ){ */

/* 	////////////////////////////////////////////////////// */
/* 	// Propose a pair of merge/split collective changes // */
/* 	////////////////////////////////////////////////////// */

/* 	// Merge ///////////////////////////////////////////// */
/* 	target = floor(gsl_rng_uniform(gen) * nnod); */
/* 	g1 = nlist[target]->inGroup; */

/* 	if(glist[g1]->size < nnod){ */

/* 	  do{ */
/* 	    target = floor(gsl_rng_uniform(gen) * nnod); */
/* 	    g2 = nlist[target]->inGroup; */
/* 	  }while( g1 == g2 ); */

/* 	  // Calculate the change of energy */
/* 	  nlink = NG2GLinksWeight(glist[g1],glist[g2]); */

/* 	  dE = 0.0; */

/* 	  dE -= (double)(2*glist[g1]->inlinksW) / (double)totallinks - */
/* 	    (double)((glist[g1]->totlinksW + glist[g1]->inlinksW) * */
/* 		     (glist[g1]->totlinksW + glist[g1]->inlinksW) ) / */
/* 	    (double)(totallinks * totallinks ); */

/* 	  dE -= (double)(2*glist[g2]->inlinksW) / (double)totallinks - */
/* 	    (double)((glist[g2]->totlinksW + glist[g2]->inlinksW) * */
/* 		     (glist[g2]->totlinksW + glist[g2]->inlinksW) ) / */
/* 	    (double)(totallinks * totallinks ); */

/* 	  dE += 2.0*(double)(glist[g1]->inlinksW +  */
/* 			     glist[g2]->inlinksW+nlink) */
/* 	    / (double)totallinks - */
/* 	    (double)(glist[g1]->totlinksW + glist[g1]->inlinksW + */
/* 		     glist[g2]->totlinksW + glist[g2]->inlinksW ) * */
/* 	    (double)(glist[g1]->totlinksW + glist[g1]->inlinksW + */
/* 		     glist[g2]->totlinksW + glist[g2]->inlinksW ) / */
/* 	    (double)(totallinks*totallinks); */

/* 	  // Accept the change according to Metroppolis */
/* 	  if( (dE >= 0.0) || ( gsl_rng_uniform(gen) < exp(dE/T) ) ){ */
/* 	    MergeGroups(glist[g1],glist[g2]); */
/* 	    energy += dE; */
/* 	  } */
/* 	} */

/* 	// Split ///////////////////////////////////////////// */
/* 	target = floor(gsl_rng_uniform(gen) * nnod); // target node */
/* 	target = nlist[target]->inGroup;    // target group */

/* 	// Look for an empty group */
/* 	g = part; */
/* 	empty = -1; */
/* 	while((g->next != NULL) && (empty < 0)){ */
/* 	  g = g->next; */
/* 	  if (g->size == 0){ */
/* 	    empty = g->label; */
/* 	  } */
/* 	} */

/* 	if (empty >= 0 ){ // if there are no empty groups, do nothing */
/* /\* 	  split = BestNetworkSplitWeight(glist[target],gen); *\/ */
/* /\* 	  split = ThermalNetworkSplitWeight(glist[target],Ti,T,gen); *\/ */
/* 	  split = ThermalPercNetworkSplitWeight(glist[target], */
/* 						Ti, T, gen); */

/* 	  // Split the group */
/* 	  nod = (split->next)->nodeList; */
/* 	  while ( nod->next != NULL ){ */
/* 	    nod = nod->next; */
/* 	    MoveNode(nlist[trans[nod->node]], */
/* 		     glist[target], */
/* 		     glist[empty]); */
/* 	  } */
/* 	  RemovePartition(split); */
/* 	  split = NULL; */

/* 	  // Try to re-merge the two groups */
/* 	  // Calculate the change of energy */
/* 	  nlink = NG2GLinksWeight(glist[target],glist[empty]); */

/* 	  dE = 0.0; */

/* 	  dE -= (double)(2*glist[target]->inlinksW) / */
/* 	    (double)totallinks - */
/* 	    (double)((glist[target]->totlinksW + */
/* 		      glist[target]->inlinksW) * */
/* 		     (glist[target]->totlinksW + */
/* 		      glist[target]->inlinksW) ) / */
/* 	    (double)(totallinks * totallinks ); */

/* 	  dE -= (double)(2*glist[empty]->inlinksW) / */
/* 	    (double)totallinks - */
/* 	    (double)((glist[empty]->totlinksW + */
/* 		      glist[empty]->inlinksW) * */
/* 		     (glist[empty]->totlinksW + */
/* 		      glist[empty]->inlinksW) ) / */
/* 	    (double)(totallinks * totallinks ); */

/* 	  dE += 2.0*(double)(glist[target]->inlinksW + */
/* 			     glist[empty]->inlinksW+nlink) */
/* 	    / (double)totallinks - */
/* 	    (double)(glist[target]->totlinksW + */
/* 		     glist[target]->inlinksW + */
/* 		     glist[empty]->totlinksW + */
/* 		     glist[empty]->inlinksW) * */
/* 	    (double)(glist[target]->totlinksW + */
/* 		     glist[target]->inlinksW + */
/* 		     glist[empty]->totlinksW + */
/* 		     glist[empty]->inlinksW) / */
/* 	    (double)(totallinks*totallinks); */

/* 	  // Accept the change according to "inverse" Metroppolis. */
/* 	  // Inverse means that the algor is applied to the split */
/* 	  // and NOT to the merge! */
/* 	  if( (dE >= 0.0) && ( gsl_rng_uniform(gen) > exp(-dE/T) ) ){ */
/* 	    MergeGroups(glist[target],glist[empty]); */
/* 	  } */
/* 	  else{ */
/* 	    energy -= dE; */
/* 	  } */

/* 	} // End of if empty */

/*       } // End of cicle2 loop */
/*     } // End of if merge == 1 */

/*     for ( i=0; i < cicle1; i++ ){ */

/*       /////////////////////////////// */
/*       // Propose an individual change */
/*       /////////////////////////////// */
/*       target = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(gsl_rng_uniform(gen) * (double)nnod); */
/*       }while(newg == oldg); */

/*       // Calculate the change of energy */
/*       inold = NodeDegreeWeightInGroup(nlist[target],glist[oldg]); */
/*       innew = NodeDegreeWeightInGroup(nlist[target],glist[newg]); */
/*       nlink = NodeDegreeWeight(nlist[target]); */

/*       dE = 0.0; */

/*       dE -= (double)(2 * glist[oldg]->inlinksW) / */
/* 	(double)totallinks -  */
/* 	(double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) * */
/* 	(double)(glist[oldg]->totlinksW+glist[oldg]->inlinksW) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE -= (double)(2 * glist[newg]->inlinksW) / */
/* 	(double)totallinks -  */
/* 	(double)(glist[newg]->totlinksW+glist[newg]->inlinksW) * */
/* 	(double)(glist[newg]->totlinksW+glist[newg]->inlinksW) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE += (double)(2*glist[oldg]->inlinksW - 2*inold) / */
/* 	(double)totallinks - */
/* 	(double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW - */
/* 		 nlink ) * */
/* 	(double)(glist[oldg]->totlinksW + glist[oldg]->inlinksW - */
/* 		 nlink ) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE += (double)(2*glist[newg]->inlinksW + 2*innew) / */
/* 	(double)totallinks - */
/* 	(double)(glist[newg]->totlinksW + glist[newg]->inlinksW + */
/* 		 nlink ) * */
/* 	(double)(glist[newg]->totlinksW + glist[newg]->inlinksW + */
/* 		 nlink ) / */
/* 	((double)totallinks * (double)totallinks); */

/*       // Accept the change according to Metroppolis */
/*       if( (dE >= 0.0) || ( gsl_rng_uniform(gen) < exp(dE/T) ) ){ */
/* 	MoveNode(nlist[target],glist[oldg],glist[newg]); */
/* 	energy += dE; */
/*       } */
/*     } */

/*     T = T * Ts; */
/*   } */

/* /\*   printf("energy = %lf\n",energy); *\/ */

/*   return CompressPart(part); */
/* } */




/* struct group *AllLevelsCommunities(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int merge, gsl_rng *gen) */
/* { */
/*   struct group *hier = NULL; */
/*   struct group *plist[2*max_size+1]; */
/*   struct group *g; */
/*   struct node_gra *nettemp = NULL; */
/*   int lastpart,npart,new_npart; */
/*   int i; */
/*   double Ti_loc; */
/*   int nnod; */

/*   for (i=0; i<(2*max_size+1); i++){ */
/*     plist[i] = NULL; */
/*   } */

/*   nnod = CountNodes(net); */
/*   hier = SACommunityIdent(net,Ti,Tf,Ts,fac,merge,gen); */
/*   plist[0] = hier; */
/*   lastpart = 0; */
/*   npart = 1; */

/*   while(npart > 0){ */

/*     new_npart = 0; */

/*     for(i=lastpart; i<(lastpart+npart); i++){ */

/*       g = plist[i]; */
/*       printf("Considering partition:\n"); */
/*       PrintGroups(plist[i]); */
/*       printf("\n"); */

/*       while(g->next != NULL){ */

/* 	g = g->next; */

/* 	// Build the network with only the nodes in the group */
/* 	nettemp = BuildNetFromGroup(g); */

/* 	if(CountNodes(nettemp) > 1){ */

/* 	  // Find the communities inside this group */
/* 	  plist[lastpart+new_npart+npart] =  */
/* 	    SACommunityIdent(nettemp,Ti*(double)nnod / */
/* 			     (double)CountNodes(nettemp), */
/* 			     Tf,Ts,fac,merge,gen); */

/* 	  // Check if the network has been split or not */
/* /\* 	  if(plist[lastpart+new_npart+npart]->next->size < g->size){ *\/ */
/* 	  if(Modularity(plist[lastpart+new_npart+npart]) > 0.00){ */
/* 	    // The network has been split! */

/* 	    printf("Modularity = %lf --- SPLITTING!\n\n",Modularity(plist[lastpart+new_npart+npart])); */

/* 	    // The pointers of the node_lis in the groups are pointing */
/* 	    // to nettemp, and they need to point to the original net! */
/* 	    MapPartToNet(plist[lastpart+new_npart+npart],net); */

/* 	    // Set this part as the offspring of the group */
/* 	    g->offspr = plist[lastpart+new_npart+npart]; */

/* 	    new_npart++; */
/* 	  } */

/* 	  else{ // Remove this partition */

/* 	    printf("Modularity = %lf --- DISREGARDING!\n\n",Modularity(plist[lastpart+new_npart+npart])); */

/* 	    RemovePartition(plist[lastpart+new_npart+npart]); */
/* 	    plist[lastpart+new_npart+npart] = NULL; */
/* 	  } */
/* 	} */

/* 	// Remove the temporal network */
/* 	RemoveGraph(nettemp); */
/* 	nettemp = NULL; */
/*       } */
/*     } */

/*     lastpart = lastpart + npart; */
/*     npart = new_npart; */
/*   } */

/*   return hier; */
/* } */
