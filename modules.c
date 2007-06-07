#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include "prng.h"

#include "tools.h"
#include "graph.h"
#include "modules.h"

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Group creation and memory allocation
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Create an empty group to be used as the partition header
// ---------------------------------------------------------------------
struct group *CreateHeaderGroup()
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

// ---------------------------------------------------------------------
// Create a group at the end of the list starting at part
// ---------------------------------------------------------------------
struct group *CreateGroup(struct group *part, int label)
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


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Partition creation
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Build a partition from a file. The file should contain
// ---------------------------------------------------------------------
struct group *FBuildPartition(FILE *inF)
{
  char label[MAX_LABEL_LENGTH];
  char *separator = "///";
  struct group *g = NULL;
  struct group *part = NULL;
  int npart = 0;

  // Create the header of the partition
  part = CreateHeaderGroup();

  // Read the file
  while (!feof(inF)) {
    g = CreateGroup(part, npart);
    npart++;
    fscanf(inF, "%s\n", &label[0]);
    while (strcmp(label, separator) != 0) {
      AddNodeToGroupSoft(g, label);
      fscanf(inF, "%s\n", &label[0]);
    }
  }

  return part;
}

// ---------------------------------------------------------------------
// For an arbitrary network, create a partition with groups of nodes
// of a given size. The nodes are placed in groups in order.
// ---------------------------------------------------------------------
struct group *CreateEquiNPartition(struct node_gra *net, int gsize)
{
  int i,j;
  int ngroups;
  struct group *g = NULL;
  struct group *part = NULL;
  struct node_gra *p = NULL;

  // Initialize stuff
  p = net;
  part = CreateHeaderGroup();
  ngroups = CountNodes(net) / gsize;

  // Create the groups and assign the nodes
  for (i=0; i<=ngroups; i++) {
    g = CreateGroup(part,i);
    for (j=0; j<gsize; j++) {
      if ((p = p->next) != NULL) {
	AddNodeToGroup(g, p);
      }
    }
  }

  // Done
  return CompressPart(part);
}

// ---------------------------------------------------------------------
// Create a partition with a given number of groups of a given
// size. The partition is not mapped to any network, and labels are
// set to integer numbers from 1 to N.
// ---------------------------------------------------------------------
struct group *CreateEquiNPartitionSoft(int ngroups, int gsize)
{
  int i, j;
  struct group *g = NULL;
  struct group *part = NULL;
  char label[MAX_LABEL_LENGTH];

  // Initialize
  part = CreateHeaderGroup();

  // Create the groups and assign the nodes
  for (i=0; i<ngroups; i++) {
    g = CreateGroup(part, i);
    for (j=0; j<gsize; j++) {
      sprintf(&label[0], "%d", (j + i*gsize + 1));
      AddNodeToGroupSoft(g, label);
    }
  }
  
  // Done
  return part;
}


/* struct group *CreatePartitionFromInGroup(struct node_gra *net) */
/* { */
/*   int i; */
/*   struct group *glist[max_size]; */
/*   struct group *part = NULL; */
/*   struct node_gra *p = NULL; */

/*   for(i=0; i<max_size; i++) */
/*     glist[i] = NULL; */

/*   part = CreateHeaderGroup(); */

/*   // Assign the nodes creating new groups when necessary */
/*   p = net; */
/*   while ((p = p->next) != NULL) { */
/*     if (glist[p->inGroup] == NULL) */
/*       glist[p->inGroup] = CreateGroup(part, p->inGroup); */
    
/*     AddNodeToGroup(glist[p->inGroup], p); */
/*   } */

/*   return part; */
/* } */


/* struct group *CreateCustomPartition(struct node_gra *net, int group[]) */
/* { */
/*   int i; */
/*   struct group *glist[max_size]; */
/*   struct group *part = NULL; */
/*   struct node_gra *p = NULL; */

/*   for(i=0; i<max_size; i++) */
/*     glist[i] = NULL; */

/*   part = CreateHeaderGroup(); */

/*   // Assign the nodes creating new groups when necessary */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */

/*     if(glist[group[p->num]] == NULL) */
/*       glist[group[p->num]] = CreateGroup(part,group[p->num]); */

/*     AddNodeToGroup(glist[group[p->num]],p); */
/*   } */

/*   return part; */
/* } */


/* struct group *CreateCustomPartitionSoft(int S,int group[]) */
/* { */
/*   int i; */
/*   struct group *glist[max_size]; */
/*   struct group *part = NULL; */

/*   for(i=0; i<max_size; i++) */
/*     glist[i] = NULL; */

/*   part = CreateHeaderGroup(); */

/*   // Assign the nodes creating new groups when necessary */
/*   for(i=0; i<S; i++){ */

/*     if(glist[group[i]] == NULL) */
/*       glist[group[i]] = CreateGroup(part,group[i]); */

/*     AddNodeToGroupSoft(glist[group[i]],i); */
/*   } */

/*   return part; */
/* } */


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Partition removal
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Free the memory allocated to a partition
// ---------------------------------------------------------------------
void RemovePartition(struct group *part)
{
  if (part->next != NULL)
    RemovePartition(part->next);

  if (part->nodeList != NULL)
    FreeAdjacencyList(part->nodeList);

  if (part->offspr != NULL)
    RemovePartition(part->offspr);

  free(part);
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Node-group functions
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Add a node to a group, pointing all pointers and updating group and
// node properties.
// ---------------------------------------------------------------------
struct node_lis *AddNodeToGroup(struct group *g, struct node_gra *node)
{
  struct node_lis *p = g->nodeList;
  int totlink, inlink;
  double totweight, inweight;

  // Go to the end of the list of nodes in the group
  while (p->next != NULL)
    p = p->next;

  // Create the node_lis and point it to the node
  p->next = (struct node_lis *)calloc(1, sizeof(struct node_lis));
  (p->next)->node = node->num;
  (p->next)->nodeLabel = (char *) calloc(MAX_LABEL_LENGTH, sizeof(char));
  strcpy((p->next)->nodeLabel, node->label);
  (p->next)->status = 0;
  (p->next)->next = NULL;
  (p->next)->ref = node;
  
  (p->next)->btw=0.0;  //initalising rush variable to 0.0;
  (p->next)->weight=0.0;  //initalising rush variable to 0.0;

  // Update the properties of the group
  g->size++;
  totlink = CountLinks(node);
  inlink = NLinksToGroup(node, g);
  totweight = NodeStrength(node);
  inweight = StrengthToGroup(node, g);
  g->totlinks += totlink - inlink;
  g->inlinks += inlink;
  g->outlinks = g->totlinks - g->inlinks;
  g->totlinksW += totweight - inweight;
  g->inlinksW += inweight;
  g->outlinksW = g->totlinksW - g->inlinksW;

  // Update the properties of the node
  node->inGroup = g->label;

  // Done
  return p->next;
}

// ---------------------------------------------------------------------
// Softly add a node to a group. Pointers do NOT point anywhere and
// there is NO updating group or node properties.
// ---------------------------------------------------------------------
struct node_lis *AddNodeToGroupSoft(struct group *g, char *label)
{
  struct node_lis *p = g->nodeList;

  // Go to the end of the list of nodes in the group
  while (p->next != NULL)
    p = p->next;

  // Create the node_lis
  p->next = (struct node_lis *)calloc(1, sizeof(struct node_lis));
  (p->next)->node = -1;
  (p->next)->nodeLabel = (char *) calloc(MAX_LABEL_LENGTH, sizeof(char));
  strcpy((p->next)->nodeLabel, label);
  (p->next)->status = 0;
  (p->next)->next = NULL;
  (p->next)->ref = NULL;

  (p->next)->btw=0.0;  //initalising rush variable to 0.0;
  (p->next)->weight=0.0;  //initalising rush variable to 0.0;

  // Actualize the properties of the group
  g->size++;

  // Done
  return p->next;
}

// ---------------------------------------------------------------------
// Remove a node from a group. Returns 1 if the node has been
// successfully removed and 0 if the node is not found in the group.
// ---------------------------------------------------------------------
int RemoveNodeFromGroup(struct group *g, struct node_gra *node)
{
  struct node_lis *p = g->nodeList;
  struct node_lis *temp;
  int totlink,inlink;
  double totweight,inweight;

  // Find the node
  while ((p->next != NULL) && ((p->next)->ref != node))
    p = p->next;
  
  if (p->next == NULL)
    return 0;
  else{
    temp = p->next;
    p->next = (p->next)->next;
    free(temp->nodeLabel);
    free(temp);

    // Update the properties of the group
    g->size--;
    totlink = CountLinks(node);
    inlink = NLinksToGroup(node, g);
    totweight = NodeStrength(node);
    inweight = StrengthToGroup(node, g);
    g->totlinks -= totlink - inlink;
    g->inlinks -= inlink;
    g->outlinks = g->totlinks - g->inlinks;
    g->totlinksW -= totweight - inweight;
    g->inlinksW -= inweight;
    g->outlinksW = g->totlinksW - g->inlinksW;

    // Done
    return 1;
  }
}

// ---------------------------------------------------------------------
// Move a node from old group to new group. Returns 1 if successful,
// and 0 if node is not in the old group.
// ---------------------------------------------------------------------
int MoveNode(struct node_gra *node, struct group *old, struct group *new)
{
  if (RemoveNodeFromGroup(old,node) == 0)
    return 0;

  AddNodeToGroup(new, node);
  return 1;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Group and partition operations
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Removes all empty groups from a partition
// ---------------------------------------------------------------------
struct group *CompressPart(struct group *part)
{
  struct group *g1;
  struct group *gtemp;

  g1 = part;
  while(g1->next != NULL){

    if(g1->next->size == 0){
      free(g1->next->nodeList);
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

// ---------------------------------------------------------------------
// Given a partition and a group label, return the group with the label
// ---------------------------------------------------------------------
struct group *GetGroup(struct group *part, int label)
{
  struct group *g=part;

  while ((g = g->next) != NULL)
    if (g->label == label)
      return g;

  return NULL;
}

// ---------------------------------------------------------------------
// Count the number of groups (after a given group, usually the header
// of the partition)
// ---------------------------------------------------------------------
int NGroups(struct group *part)
{
  int ngroup = 0;

  while ((part = part->next) != NULL)
    ngroup++;

  return ngroup;
}

// ---------------------------------------------------------------------
// Count the number of non-empty groups (after a given group, usually
// the header of the partition)
// ---------------------------------------------------------------------
int NNonEmptyGroups(struct group *part)
{
  int ngroup = 0;

  while ((part = part->next) != NULL)
    if (part->size > 0)
      ngroup++;

  return ngroup;
}

// ---------------------------------------------------------------------
// Count the number of nodes in all groups (after a given group,
// usually the header of the partition)
// ---------------------------------------------------------------------
int PartitionSize(struct group *part)
{
  int nnod = 0;

  while ((part = part->next) != NULL)
    nnod += part->size;

  return nnod;
}

// ---------------------------------------------------------------------
// Removes all links to and from the nodes in a group
// ---------------------------------------------------------------------
void RemoveWithinGroupLinks(struct group *g, int symmetric_sw)
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

// ---------------------------------------------------------------------
// Removes all links between groups
// ---------------------------------------------------------------------
void RemoveBetweenGroupLinks(struct group *part, int symmetric_sw)
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

// ---------------------------------------------------------------------
// 
// ---------------------------------------------------------------------
/* double **BlockModel(FILE *outf, */
/* 		    struct group *part, */
/* 		    char type_sw, */
/* 		    int list_sw) */
/* { */
/*   struct group *g1, *g2; */
/*   int links = 0; */
/*   int nodes = 0; */
/*   double prob; */
/*   double bij, av, sig; */

/*   // Count the total number of nodes and links, and the average */
/*   // linking probability */
/*   g1 = part; */
/*   while ((g1 = g1->next) != NULL) { */
/*     if (g1->size > 0) { */
/*       nodes += g1->size; */
/*       links += g1->inlinks; */
/*       g2 = g1; */
/*       while ((g2 = g2->next) != NULL) { */
/* 	if (g2->size > 0) { */
/* 	  links += CountG2GLinks(g1, g2); */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   prob = (double)links / (double)(nodes * (nodes - 1)); */

/*   // Calculate and print the blockmodel */
/*   g1 = part; */
/*   while ((g1 = g1->next) != NULL) { */
/*     if (g1->size > 0) { */
/*       g2 = part; */
/*       while ((g2 = g2->next) != NULL) { */
/* 	if (g2->size > 0) { */

/* 	  // Calculate the matrix element */
/* 	  switch (type_sw) { */
/* 	  case 'n': */
/* 	    bij = (double)CountG2GLinks(g1, g2); */
/* 	    break; */
/* 	  case 'f': */
/* 	    bij = (double)CountG2GLinks(g1, g2)/ (double)links; */
/* 	    break; */
/* 	  case 'p': */
/* 	    if (g1 == g2)  */
/* 	      bij = (double)CountG2GLinks(g1, g2) /  */
/* 		(double)(g1->size * (g1->size - 1)); */
/* 	    else */
/* 	      bij = (double)CountG2GLinks(g1, g2) /  */
/* 		(double)(g1->size * g2->size); */
/* 	    break; */
/* 	  case 'e': */
/* 	    bij = (double)CountG2GLinks(g1, g2) / (double)links - */
/* 	      (double)((g1->inlinks + g1->totlinks) * */
/* 		       (g2->inlinks + g2->totlinks)) /  */
/* 	      ((double)(4 * links * links)); */
/* 	    break; */
/* 	  case 'z': */
/* 	    if (g1 == g2) { */
/* 	      av = prob * (double)(g1->size * (g1->size - 1)); */
/* 	      sig = sqrt((double)(g1->size * (g1->size - 1)) *  */
/* 			 prob * (1.0 - prob)); */
/* 	      bij = ((double)CountG2GLinks(g1, g2) - av) / sig; */
/* 	    } */
/* 	    else { */
/* 	      av = prob * (double)(g1->size * g2->size); */
/* 	      sig = sqrt((double)(g1->size * g2->size) * prob * (1.0 - prob)); */
/* 	      bij = ((double)CountG2GLinks(g1, g2) - av) / sig; */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   return; */
/* } */

// ---------------------------------------------------------------------
// Count the number of links from a node to a given group
// ---------------------------------------------------------------------
int NLinksToGroup(struct node_gra* node, struct group *g)
{
  struct node_lis *nei = node->neig;
  int inlink = 0;

  while ((nei = nei->next) != NULL)
    if ((nei->ref)->inGroup == g->label)
      inlink++;

  return inlink;
}

// ---------------------------------------------------------------------
// Weight of links from a node to a given group
// ---------------------------------------------------------------------
double StrengthToGroup(struct node_gra* node, struct group *g)
{
  struct node_lis *nei = node->neig;
  double inlink = 0.0;

  while ((nei = nei->next) != NULL)
    if ((nei->ref)->inGroup == g->label)
      inlink += nei->weight;
  
  return inlink;
}

// ---------------------------------------------------------------------
// Count the number of links between a pair of groups
// ---------------------------------------------------------------------
int NG2GLinks(struct group *g1, struct group *g2)
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

// ---------------------------------------------------------------------
// Total weight of links between a pair of groups
// ---------------------------------------------------------------------
double NG2GLinksWeight(struct group *g1, struct group *g2)
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

// ---------------------------------------------------------------------
// Move all nodes in group g1 to group g2
// ---------------------------------------------------------------------
void MergeGroups(struct group *g1, struct group *g2)
{
  struct node_lis *p = g1->nodeList;

  while (p->next != NULL)
    MoveNode((p->next)->ref, g1, g2);
  return;
}

// ---------------------------------------------------------------------
// Creates a copy of a group and puts it at the end of the list of
// groups hanging from copy_root
// ---------------------------------------------------------------------
struct group *CopyGroup(struct group *copy_root, struct group *g)
{
  struct group *copy;
  struct node_lis *node;

  copy = CreateGroup(copy_root, g->label);

  // Copy the trivial fields
  copy->coorX = g->coorX;
  copy->coorY = g->coorY;
  copy->coorZ = g->coorZ;

  // Copy the node list
  node = g->nodeList;
  while ( node->next != NULL ) {
    node = node->next;
    AddNodeToGroup(copy, node->ref);
  }

  // Some recursivity is in order here :(
  // AND NEEDS TO BE TESTED BEFORE USING!!!!!!!!!!
  if ( g->offspr != NULL )
    copy->offspr = CopyPartition(g->offspr);

  return copy;
}

// ---------------------------------------------------------------------
// Creates a copy of a whole partition. CAUTION!!!! The network MUST
// be mapped to the copy or the original partition if they are to be
// used after the copy has been created!!!!
// ---------------------------------------------------------------------
struct group *CopyPartition(struct group *original)
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

// ---------------------------------------------------------------------
// Creates a network that contains only the nodes in one of the groups
// of a partition
// ---------------------------------------------------------------------
struct node_gra *BuildNetFromPart(struct group *part)
{
  struct node_gra *net = NULL;
  struct node_gra *last = NULL;
  struct node_gra *new = NULL;
  struct node_lis *p = (part->nodeList);

  net = CreateHeaderGraph();
  last = net;

  while ((p = p->next) != NULL) {
    new = CreateNodeGraph(last, p->ref->label);
    new->state = p->ref->state;
    new->coorX = p->ref->coorX;
    new->coorY = p->ref->coorY;
    new->coorZ = p->ref->coorZ;
    new->ivar1 = p->ref->ivar1;
    new->dvar1 = p->ref->dvar1;
    new->inGroup = p->ref->inGroup;
    last = new;
    CopyAdjacencyList(p->ref, new);
  }
  
  CleanAdjacencies(net);
  RewireAdjacency(net);

  return net;
}


// ---------------------------------------------------------------------
// Creates a network that contains only the nodes in one of the groups
// of a partition AND its neighbors
// ---------------------------------------------------------------------
struct node_gra *BuildNetFromPartNeig(struct group *part)
{
  struct node_gra *net = NULL;
  struct node_gra *last = NULL;
  struct node_gra *new = NULL;
  struct node_lis *p = part->nodeList;
  struct node_lis *nei = NULL;

  net = CreateHeaderGraph();
  last = net;

  // Loop over nodes in the partition
  while ((p = p->next) != NULL){
    // Add the node
    new = CreateNodeGraph(last, p->ref->label);
    new->state = p->ref->state;
    new->coorX = p->ref->coorX;
    new->coorY = p->ref->coorY;
    new->coorZ = p->ref->coorZ;
    new->inGroup = p->ref->inGroup;
    last = new;
    CopyAdjacencyList(p->ref, new);

    // Add the neighbors of the node
    nei = p->ref->neig;
    while ((nei = nei->next) != NULL) {
      // If the node does not exist and is in another module, add it
      // to the network
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
  RewireAdjacency(net);

  return net;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Network-partition operations
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Reset the inGroup attribute of all the nodes in a network 
// ---------------------------------------------------------------------
void ResetNetGroup(struct node_gra *net)
{
  while ((net = net->next) != NULL)
    net->inGroup = -1;

  return;
}

// ---------------------------------------------------------------------
// Given a partition and a network, set all the pointers and all the
// group and node attributes to the right values.
// ---------------------------------------------------------------------
void MapPartToNet(struct group *part, struct node_gra *net)
{
  struct group *g = NULL;
  struct node_lis *nod;
  int totlink, inlink;
  double totlinkW, inlinkW;
  void *nodeDict = NULL;
  struct node_gra *node = NULL;

  // Reset the group of all nodes
  ResetNetGroup(net);

  // Create the nodeDict for fast access to nodes by label
  nodeDict = MakeLabelDict(net);

  // Go through the groups, reset attributes, and point pointers
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
      // Get the node_gra by label
      node = GetNodeDict(nod->nodeLabel, nodeDict);

      // Update the properties of the group
      nod->ref = node;
      nod->node = node->num;
      // Update the properties of the node
      node->inGroup = g->label;
    }
  }

  // Count links in groups
  g = part;
  while ((g = g->next) != NULL) {
    nod = g->nodeList;
    while ((nod = nod->next) != NULL) {
      totlink = CountLinks(nod->ref);
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

  // Free memory allocated locally
  FreeLabelDict(nodeDict);

  // Done
  return;
}

// ---------------------------------------------------------------------
// Same as MapPartToNet but ONLY the network (not the partition) is
// updated. This enables one to map a partition to networks that are a
// subset of the original network on which the partition is based.
// ---------------------------------------------------------------------
void MapPartToNetSoft(struct group *part, struct node_gra *net)
{
  struct group *g = NULL;
  struct node_lis *nod;
  void *nodeDict = NULL;
  struct node_tree *treeNode = NULL;
  struct node_gra *node = NULL;

  // Reset the group of all nodes
  ResetNetGroup(net);

  // Create the nodeDict for fast access to nodes by label
  nodeDict = MakeLabelDict(net);

  // Go through the groups, reset attributes, and point pointers
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
      // Get the node_gra by label
      treeNode = CreateNodeTree();
      strcpy(treeNode->label, nod->nodeLabel);
      node = (*(struct node_tree **)tfind((void *)treeNode,
					  &nodeDict,
					  NodeTreeLabelCompare))->ref;
      FreeNodeTree(treeNode, preorder, 0);
      // Update the properties of the node
      node->inGroup = g->label;
    }
  }

  return;
}

// ---------------------------------------------------------------------
// Given an undirected network, create a partition in which each group
// corresponds to an isolated cluster of the network.
// ---------------------------------------------------------------------
struct group *ClustersPartition(struct node_gra *net)
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

  // Initialize some variables
  part = CreateHeaderGroup();
  size = 0;
  ResetNetGroup(net);
  nnod = CountNodes(net);
  list = CreateHeaderList();

  // Start the search for clusters
  last = net;
  do{
    // Find the first unclassified node...
    p = last->next;
    while (p->inGroup >= 0)
      p = p->next;
    last = p;

    // Create a new group in the partition and add the node
    thisgroup = CreateGroup(part, groupcoun++);
    AddNodeToGroup(thisgroup, p);

    // ...and enqueue it
    ResetNodesState(net);
    d = 0;
    Enqueue(p, list, list, &size, d);
    anod++;
    lp = list;

    // Update the list with successive neighbors
    do{
      size_ant = size;
      lp = RenewQueue(list, lp, &size, d);
      lp2 = lp;

      while (lp2->next != NULL) {
	lp2 = lp2->next;
	anod++;

	// Add the node to the group
	AddNodeToGroup(thisgroup, lp2->ref);
      }
      d++;
    } while(size != size_ant);

    ClearList(list, &size);
  } while(anod < nnod);

  // Free memory
  free(list);

  // Done
  return part;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Group and partition output
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Print the groups in a partition. If list_sw == 0, each group is
// printed in a row with some group info (nodes, links,
// etc.). Otherwise, each node is printed in a row and groups are
// separated by a line with '//'
// ---------------------------------------------------------------------
void FPrintPartition(FILE *outf, struct group *partition, int list_sw)
{
  struct group *g = partition;
  struct node_lis *p;

  while ((g = g->next) != NULL) {
    if(g->size > 0){

      // Print group info for the verbose (non-list) mode
      if (list_sw == 0) {
	fprintf(outf, "%d %d %d %d %d %lf %lf %lf ---",
		g->label+1, g->size, g->totlinks, g->inlinks, g->outlinks,
		g->totlinksW, g->inlinksW, g->outlinksW);
      }
      
      // Print nodes
      p = g->nodeList;
      while ((p = p->next) != NULL) {
	if (list_sw == 0)
	  fprintf(outf, " %s", p->nodeLabel);
	else
	  fprintf(outf, "%s\n", p->nodeLabel);
      }

      // End group
      if (list_sw == 0)
	fprintf(outf, "\n");
      else
	fprintf(outf, "//\n");
    }
  } // End of loop over groups in the partition

  return;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Module indentification
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Calculate the modularity of a partition
// ---------------------------------------------------------------------
double Modularity(struct group *part)
{
  struct group *g = part;
  int links2 = 0;
  double modul = 0.0;

  // Calculate the number of links (times 2)
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

// ---------------------------------------------------------------------
// Calculate the weighted modularity of a partition
// ---------------------------------------------------------------------
double ModularityWeight(struct group *part)
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

// ---------------------------------------------------------------------
// Given a group, SAGroupSplit returns a split of the nodes into two
// subgroups. It uses SA, so initial and final temperature must be
// provided. If cluster_sw == 1, then the function checks first
// whether the group is disconnected; if it is, then it returns (with
// some probability) a partition along the lines of the existing
// clusters without any SA (with complementary probability, it still
// does the SA).
// ---------------------------------------------------------------------
struct group *SAGroupSplit(struct group *targ,
			   double Ti, double Tf, double Ts,
			   int cluster_sw,
			   struct prng *gen)
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

  nodeList = (struct node_gra**)malloc(targ->size * sizeof(struct node_gra *));
  glist[0] = NULL;
  glist[1] = NULL;

  // Build a network from the nodes in the target group
  net = BuildNetFromPart(targ);

  // Check if the network is connected
  split = ClustersPartition(net);
  ngroups = NGroups(split);

  if (
      cluster_sw == 1 &&         // cluster switch is on
      ngroups > 1 &&             // network is not connected
      prng_get_next(gen) < prob  // with some probability
      ) {
    
    // Merge groups randomly until only two are left
    while (ngroups > 2) {
      // Select two random groups
      g1 = ceil(prng_get_next(gen)* (double)ngroups);
      do {
	g2 = ceil(prng_get_next(gen)* (double)ngroups);
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

  else { // Network is connected or we want to ignore the clusters
    // Remove SCS partition
    RemovePartition(split);
    ResetNetGroup(net);

    // Create the groups
    split = CreateHeaderGroup();
    glist[0] = CreateGroup(split, 0);
    glist[1] = CreateGroup(split, 1);

    // Randomly assign the nodes to the groups
    p = net;
    while ((p = p->next) != NULL) {
      nodeList[nnod] = p;
      totallinks += CountLinks(p);
      nnod++;
      
      des = floor(prng_get_next(gen) * 2.0);
      AddNodeToGroup(glist[des], p);
    }
    totallinks /= 2;
    
    // Do the SA to "optimize" the splitting
    if (totallinks > 0) {
      T = Ti;
      while (T > Tf) {
	
	for (i=0; i<nnod; i++) {
	  target = floor(prng_get_next(gen) * (double)nnod);
	  oldg = nodeList[target]->inGroup;
	  if(oldg == 0)
	    newg = 1;
	  else
	    newg = 0;
	
	  // Calculate the change of energy
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

	  // Accept the move according to Metropolis
	  if( (dE >=0.0) || (prng_get_next(gen) < exp(dE/T)) ){
	    MoveNode(nodeList[target],glist[oldg],glist[newg]);
	    energy += dE;
	  }
	}

	T = T * Ts;
      } // End of temperature loop
    } // End if totallinks > 0
  }

  // Free memory
  RemoveGraph(net);
  free(nodeList);

  // Done
  return split;
}


/* struct group *ThermalNetworkSplitWeight(struct group *targ,double Ti, double Tf, struct prng *gen) */
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
/*   net = BuildNetFromPart(targ); */

/*   // Create the groups */
/*   split = CreateHeaderGroup(); */
/*   glist[0] = CreateGroup(split,0); */
/*   glist[1] = CreateGroup(split,1); */

/*   // Randomly assign the nodes to the groups */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     nlist[nnod] = p; */
/*     totallinks += CountLinksWeight(p); */
/*     nnod++; */

/*     des = floor(prng_get_next(gen)*2.0); */
/*     AddNodeToGroup(glist[des],p); */
/*   } */

/*   totallinks /= 2.0; */

/*   // Do the SA to "optimize" the splitting */
/*   T = Ti; */
/*   while( T > Tf){ */

/*     for (i=0; i< nnod; i++){ */
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       if(oldg == 0) */
/* 	newg = 1; */
/*       else */
/* 	newg = 0; */

/*       // Calculate the change of energy */
/*       inold = CountLinksWeightInGroup(nlist[target],glist[oldg]); */
/*       innew = CountLinksWeightInGroup(nlist[target],glist[newg]); */
/*       nlink = CountLinksWeight(nlist[target]); */

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
/*       if( (dE >= 0.0) || (prng_get_next(gen) < exp(dE/T)) ){ */
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
/* 					    struct prng *gen) */
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
/*   net = BuildNetFromPart(targ); */

/*   // Check if the network is connected */
/*   split = ClustersPartition(net); */
/*   ngroups = CountGroups(split); */

/*   if ( ngroups > 1 && prng_get_next(gen) < prob) { // Network is not */
/* 						   // connected */
    
/*     // Merge groups randomly until only two are left */
/*     while (ngroups > 2) { */
/*       // Select two random groups */
/*       g1 = ceil(prng_get_next(gen)* (double)ngroups); */
/*       do { */
/* 	g2 = ceil(prng_get_next(gen)* (double)ngroups); */
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
/*       totallinks += CountLinksWeight(p); */
/*       nnod++; */
      
/*       des = floor(prng_get_next(gen)*2.0); */
/*       AddNodeToGroup(glist[des],p); */
/*     } */

/*     totallinks /= 2.0; */

/*     // Do the SA to "optimize" the splitting */
/*     if ( totallinks > 0 ) { */
/*       T = Ti; */
/*       while( T > Tf){ */
      
/* 	for (i=0; i< nnod; i++){ */
/* 	  target = floor(prng_get_next(gen) * (double)nnod); */
/* 	  oldg = nlist[target]->inGroup; */
/* 	  if(oldg == 0) */
/* 	    newg = 1; */
/* 	  else */
/* 	    newg = 0; */

/* 	  // Calculate the change of energy */
/* 	  inold = CountLinksWeightInGroup(nlist[target],glist[oldg]); */
/* 	  innew = CountLinksWeightInGroup(nlist[target],glist[newg]); */
/* 	  nlink = CountLinksWeight(nlist[target]); */
	  
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
/* 	  if( (dE >= 0.0) || (prng_get_next(gen) < exp(dE/T)) ){ */
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
/* struct group *SACommunityIdent(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int merge, struct prng *gen) */
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
/*   double energy = 0.0, dE = 0.0; */
/*   double T; */
/*   int g1,g2; */
/*   double energyant = 0.0; */
/*   int count = 0, limit = 25; // to stop the search if the energy does */
/*                              // not change */
/*   int cicle1,cicle2; */
/*   int trans[maxim_int]; */
/*   struct group *best_part = NULL; */
/*   double best_E = -100.0; */

/*   // Create the groups and assign each node to one group */
/*   nnod = CountNodes(net); */
/*   part = CreateHeaderGroup(); */
/*   p = net->next; */
/*   ResetNetGroup(net); // All nodes reset to group -1 */

/*   nlist[0] = p; */
/*   trans[p->num] = 0; */
/*   glist[0] = CreateGroup(part,0); */
/*   AddNodeToGroup(glist[0],p); */
/*   totallinks += CountLinks(p); */
  
/*   for( i=1; i<nnod; i++ ) { */
/*     p = p->next; */

/*     nlist[i] = p; */
/*     trans[p->num] = i; */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*     AddNodeToGroup(glist[i],p); */
/*     totallinks += CountLinks(p); */
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
/*   energy = Modularity(part); */

/*   while( T > Tf && count < limit){ */

/* /\*     PrintGroups(part); *\/ */
/* /\*     printf("%g %lf %lf %g\n",1.0/T, energy, Modularity(part), T); *\/ */
/* /\*     printf("%g %lf %g %g\n",1.0/T, energy, T, AverageGroupSize(part) ); *\/ */
/*     printf("%g %lf %g\n",1.0/T, energy, T); */

/*     for ( i=0; i < cicle1; i++ ){ */
      
/*       /////////////////////////////// */
/*       // Propose an individual change */
/*       /////////////////////////////// */
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(prng_get_next(gen) * (double)nnod); */
/*       }while(newg == oldg); */

/* /\*       printf("%d: %d->%d?\n", target+1, oldg+1, newg+1); *\/ */

/*       // Calculate the change of energy */
/*       inold = NLinksToGroup(nlist[target],glist[oldg]); */
/*       innew = NLinksToGroup(nlist[target],glist[newg]); */
/*       nlink = CountLinks(nlist[target]); */

/* /\*       printf("--2\n"); *\/ */

/*       dE = 0.0; */

/*       dE -= (double)(2 * glist[oldg]->inlinks) / */
/* 	(double)totallinks -  */
/* 	(double)(glist[oldg]->totlinks+glist[oldg]->inlinks) * */
/* 	(double)(glist[oldg]->totlinks+glist[oldg]->inlinks) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE -= (double)(2 * glist[newg]->inlinks) / */
/* 	(double)totallinks -  */
/* 	(double)(glist[newg]->totlinks+glist[newg]->inlinks) * */
/* 	(double)(glist[newg]->totlinks+glist[newg]->inlinks) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE += (double)(2*glist[oldg]->inlinks - 2*inold) / */
/* 	(double)totallinks - */
/* 	(double)(glist[oldg]->totlinks + glist[oldg]->inlinks - */
/* 		 nlink ) * */
/* 	(double)(glist[oldg]->totlinks + glist[oldg]->inlinks - */
/* 		 nlink ) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE += (double)(2*glist[newg]->inlinks + 2*innew) / */
/* 	(double)totallinks - */
/* 	(double)(glist[newg]->totlinks + glist[newg]->inlinks + */
/* 		 nlink ) * */
/* 	(double)(glist[newg]->totlinks + glist[newg]->inlinks + */
/* 		 nlink ) / */
/* 	((double)totallinks * (double)totallinks); */

/* /\*       printf("--4\n"); *\/ */

/*       // Accept the change according to Metroppolis */
/*       if( (dE >= 0.0) || ( prng_get_next(gen) < exp(dE/T) ) ){ */
/* 	MoveNode(nlist[target],glist[oldg],glist[newg]); */
/* 	energy += dE; */
/*       } */

/* /\*       printf("--5\n"); *\/ */
/*     } */

/* /\*     printf("Done\n"); *\/ */

/*     if (merge == 1) { */

/*       for ( i=0; i < cicle2; i++ ){ */

/* 	////////////////////////////////////////////////////// */
/* 	// Propose a pair of merge/split collective changes // */
/* 	////////////////////////////////////////////////////// */

/* 	// Merge ///////////////////////////////////////////// */
/* 	target = floor(prng_get_next(gen) * (double)nnod); */
/* 	g1 = nlist[target]->inGroup; */
	
/* 	if(glist[g1]->size < nnod){ */

/* 	  do{ */
/* 	    target = floor(prng_get_next(gen) * (double)nnod); */
/* 	    g2 = nlist[target]->inGroup; */
/* 	  }while( g1 == g2 ); */

/* /\* 	  printf("merge %d-%d\n"); *\/ */
	
/* 	  // Calculate the change of energy */
/* 	  nlink = NG2GLinks(glist[g1],glist[g2]); */
	  
/* 	  dE = 0.0; */
	  
/* 	  dE -= (double)(2*glist[g1]->inlinks) / (double)totallinks - */
/* 	    ((double)(glist[g1]->totlinks + glist[g1]->inlinks) * */
/* 	     (double)(glist[g1]->totlinks + glist[g1]->inlinks) ) / */
/* 	    ((double)totallinks * (double)totallinks ); */
	
/* 	  dE -= (double)(2*glist[g2]->inlinks) / (double)totallinks - */
/* 	    ((double)(glist[g2]->totlinks + glist[g2]->inlinks) * */
/* 	     (double)(glist[g2]->totlinks + glist[g2]->inlinks) ) / */
/* 	    ((double)totallinks * (double)totallinks ); */

/* 	  dE += 2.0*(double)(glist[g1]->inlinks + */
/* 			     glist[g2]->inlinks+nlink) */
/* 	    / (double)totallinks - */
/* 	    (double)(glist[g1]->totlinks + glist[g1]->inlinks + */
/* 		     glist[g2]->totlinks + glist[g2]->inlinks ) * */
/* 	    (double)(glist[g1]->totlinks + glist[g1]->inlinks + */
/* 		     glist[g2]->totlinks + glist[g2]->inlinks ) / */
/* 	    ((double)totallinks * (double)totallinks); */
	
/* 	  // Accept the change according to Metroppolis */
/* 	  if( (dE >= 0.0) || ( prng_get_next(gen) < exp(dE/T) ) ){ */
/* 	    MergeGroups(glist[g1],glist[g2]); */
/* 	    energy += dE; */
/* 	  } */
/* 	} */

/* 	// Split ///////////////////////////////////////////// */
/* 	target = floor(prng_get_next(gen)*(double)nnod); //target node */
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
/* /\* 	  split = BestNetworkSplit(glist[target],gen); *\/ */
/* /\* 	  split = ThermalNetworkSplit(glist[target],Ti,T,gen); *\/ */
/* 	  split = ThermalPercNetworkSplit(glist[target],Ti,T,gen); */

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
/* 	  nlink = NG2GLinks(glist[target],glist[empty]); */
	
/* 	  dE = 0.0; */
	
/* 	  dE -= (double)(2*glist[target]->inlinks) / */
/* 	    (double)totallinks - */
/* 	    ((double)(glist[target]->totlinks+glist[target]->inlinks) * */
/* 	     (double)(glist[target]->totlinks+glist[target]->inlinks))/ */
/* 	    ((double)totallinks * (double)totallinks ); */
	
/* 	  dE -= (double)(2*glist[empty]->inlinks) / */
/* 	    (double)totallinks - */
/* 	    ((double)(glist[empty]->totlinks + glist[empty]->inlinks) * */
/* 	     (double)(glist[empty]->totlinks + glist[empty]->inlinks))/ */
/* 	    ((double)totallinks * (double)totallinks ); */

/* 	  dE += 2.0*(double)(glist[target]->inlinks+ */
/* 			     glist[empty]->inlinks+nlink) */
/* 	    / (double)totallinks - */
/* 	    (double)(glist[target]->totlinks + glist[target]->inlinks */
/* 		     + glist[empty]->totlinks + */
/* 		     glist[empty]->inlinks ) * */
/* 	    (double)(glist[target]->totlinks + glist[target]->inlinks + */
/* 		     glist[empty]->totlinks + glist[empty]->inlinks ) / */
/* 	    ((double)totallinks * (double)totallinks); */
	
/* 	  // Accept the change according to "inverse" Metroppolis. */
/* 	  // Inverse means that the algor is applied to the split and */
/* 	  // NOT to the merge! */
/* 	  if( (dE > 0.0) && ( prng_get_next(gen) > exp(-dE/T) ) ){ */
/* 	    MergeGroups(glist[target],glist[empty]); */
/* 	  } */
/* 	  else{ */
/* 	    energy -= dE; */
/* 	  } */
	  
/* 	} // End of if empty >= 0 */

/*       } // End of cicle 2 loop */

/*     } // End of if merge==1 */


/*     // Update the no-change counter */
/*     if (energy == energyant) { */
/*       count++; */
      
/*       // If the program is ready to stop (count==limit) but the */
/*       // current partition is not the best one so far, replace the */
/*       // current partition by the best one and continue from there. */
/*       if ( (count == limit) && (energy < best_E) ) { */
/* 	printf("# Resetting partition\n"); */
/* 	RemovePartition(part); */
/* 	g = part = CopyPartition(best_part); */
	
/* 	while ( g->next != NULL ) { */
/* 	  g = g->next; */
/* 	  glist[g->label] = g; */
/* 	} */

/* 	MapPartToNet(part, net); */
/* 	energy = best_E; */
/* 	count = 0; */
/*       } */
/*     } */
/*     else{ */
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
/*       MapPartToNet(part, net); // MUST DO this after copying a part! */
/*       best_E = energy; */
/*     } */

/*     // Uptade the temperature */
/*     T = T * Ts; */

/*   } // End of simulated annealing */

/*   RemovePartition(best_part); */
/*   return CompressPart(part); */
/* } */


/* // merge = 0 => No group merging */
/* // Same as SACommunityIdent but the maximum number of modules */
/* // is specified */
/* struct group *SACommunityIdentNMod(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int merge, int nmod, struct prng *gen) */
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


/*   // Create the groups the nodes randomly */
/*   nnod = CountNodes(net); */
/*   part = CreateHeaderGroup(); */
/*   p = net->next; */
/*   ResetNetGroup(net); // All nodes reset to group -1 */

/*   // Create the groups */
/*   for (i=0; i<nmod; i++){ */
/*     glist[i] = CreateGroup(part,i); */
/*   } */

/*   // Assign the nodes randomly */
/*   nlist[0] = p; */
/*   trans[p->num] = 0; */
/*   target = floor( prng_get_next(gen) * (double)nmod ); */
/*   AddNodeToGroup(glist[target],p); */
/*   totallinks += CountLinks(p); */

/*   for( i=1; i<nnod; i++ ) { */
/*     p = p->next; */

/*     nlist[i] = p; */
/*     trans[p->num] = i; */
/*     target = floor( prng_get_next(gen) * (double)nmod ); */
/*     AddNodeToGroup(glist[target],p); */
/*     totallinks += CountLinks(p); */
/*   } */

/*   // Number of iterations at each temperature */
/*   if (fac*(double)(nnod*nmod) < 10) */
/*     cicle1 = 10; */
/*   else */
/*     cicle1 = floor(fac*(double)(nnod*nmod)); */

/*   if (fac*(double)nmod < 2) */
/*     cicle2 = 2; */
/*   else */
/*     cicle2 = floor(fac*(double)nmod); */

/*   // Do the simulated annealing */
/*   T = Ti; */
/*   energy = Modularity(part); */

/*   while( T > Tf && count < limit){ */

/*     if (energy == energyant) */
/*       count++; */
/*     else{ */
/*       energyant = energy; */
/*       count = 0; */
/*     } */

/* /\*     PrintGroups(part); *\/ */
/* /\*     printf("%g %lf %lf %g\n",1.0/T, energy, Modularity(part), T); *\/ */
/* /\*     printf("%g %lf %g\n",1.0/T, energy, T); *\/ */

/*     for ( i=0; i < cicle1; i++ ){ */
      
/*       /////////////////////////////// */
/*       // Propose an individual change */
/*       /////////////////////////////// */
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(prng_get_next(gen) * (double)nmod); */
/*       }while(newg == oldg); */

/*       // Calculate the change of energy */
/*       inold = NLinksToGroup(nlist[target],glist[oldg]); */
/*       innew = NLinksToGroup(nlist[target],glist[newg]); */
/*       nlink = CountLinks(nlist[target]); */

/*       dE = 0.0; */

/*       dE -= (double)(2 * glist[oldg]->inlinks) / */
/* 	(double)totallinks -  */
/* 	(double)(glist[oldg]->totlinks+glist[oldg]->inlinks) * */
/* 	(double)(glist[oldg]->totlinks+glist[oldg]->inlinks) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE -= (double)(2 * glist[newg]->inlinks) / */
/* 	(double)totallinks -  */
/* 	(double)(glist[newg]->totlinks+glist[newg]->inlinks) * */
/* 	(double)(glist[newg]->totlinks+glist[newg]->inlinks) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE += (double)(2*glist[oldg]->inlinks - 2*inold) / */
/* 	(double)totallinks - */
/* 	(double)(glist[oldg]->totlinks + glist[oldg]->inlinks - */
/* 		 nlink ) * */
/* 	(double)(glist[oldg]->totlinks + glist[oldg]->inlinks - */
/* 		 nlink ) / */
/* 	((double)totallinks * (double)totallinks); */

/*       dE += (double)(2*glist[newg]->inlinks + 2*innew) / */
/* 	(double)totallinks - */
/* 	(double)(glist[newg]->totlinks + glist[newg]->inlinks + */
/* 		 nlink ) * */
/* 	(double)(glist[newg]->totlinks + glist[newg]->inlinks + */
/* 		 nlink ) / */
/* 	((double)totallinks * (double)totallinks); */

/*       // Accept the change according to Metroppolis */
/*       if( (dE >= 0.0) || ( prng_get_next(gen) < exp(dE/T) ) ){ */
/* 	MoveNode(nlist[target],glist[oldg],glist[newg]); */
/* 	energy += dE; */
/*       } */
/*     } */

/*     if (merge == 1) { */

/*       for ( i=0; i < cicle2; i++ ){ */

/* 	////////////////////////////////////////////////////// */
/* 	// Propose a pair of merge/split collective changes // */
/* 	////////////////////////////////////////////////////// */

/* 	// Merge ///////////////////////////////////////////// */
/* 	target = floor(prng_get_next(gen) * (double)nmod); */
/* 	g1 = nlist[target]->inGroup; */
	
/* 	if(glist[g1]->size < nnod){ */

/* 	  do{ */
/* 	    target = floor(prng_get_next(gen) * (double)nnod); */
/* 	    g2 = nlist[target]->inGroup; */
/* 	  }while( g1 == g2 ); */
	
/* 	  // Calculate the change of energy */
/* 	  nlink = NG2GLinks(glist[g1],glist[g2]); */
	  
/* 	  dE = 0.0; */
	  
/* 	  dE -= (double)(2*glist[g1]->inlinks) / (double)totallinks - */
/* 	    (double)((glist[g1]->totlinks + glist[g1]->inlinks) * */
/* 		     (glist[g1]->totlinks + glist[g1]->inlinks) ) / */
/* 	    (double)(totallinks * totallinks ); */
	
/* 	  dE -= (double)(2*glist[g2]->inlinks) / (double)totallinks - */
/* 	    (double)((glist[g2]->totlinks + glist[g2]->inlinks) * */
/* 		     (glist[g2]->totlinks + glist[g2]->inlinks) ) / */
/* 	    (double)(totallinks * totallinks ); */

/* 	  dE += 2.0*(double)(glist[g1]->inlinks +  */
/* 			     glist[g2]->inlinks+nlink) */
/* 	    / (double)totallinks - */
/* 	    (double)(glist[g1]->totlinks + glist[g1]->inlinks + */
/* 		     glist[g2]->totlinks + glist[g2]->inlinks ) * */
/* 	    (double)(glist[g1]->totlinks + glist[g1]->inlinks + */
/* 		     glist[g2]->totlinks + glist[g2]->inlinks ) / */
/* 	    (double)(totallinks*totallinks); */
	
/* 	  // Accept the change according to Metroppolis */
/* 	  if( (dE >= 0.0) || ( prng_get_next(gen) < exp(dE/T) ) ){ */
/* 	    MergeGroups(glist[g1],glist[g2]); */
/* 	    energy += dE; */
/* 	  } */
/* 	} */

/* 	// Split ///////////////////////////////////////////// */
/* 	target = floor(prng_get_next(gen) * (double)nnod); // target node */
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
/* /\* 	  split = BestNetworkSplit(glist[target],gen); *\/ */
/* 	  split = ThermalNetworkSplit(glist[target],Ti,T,gen); */

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
/* 	  nlink = NG2GLinks(glist[target],glist[empty]); */
	
/* 	  dE = 0.0; */
	
/* 	  dE -= (double)(2*glist[target]->inlinks) / (double)totallinks - */
/* 	    (double)((glist[target]->totlinks + glist[target]->inlinks) * */
/* 		     (glist[target]->totlinks + glist[target]->inlinks) ) / */
/* 	    (double)(totallinks * totallinks ); */
	
/* 	  dE -= (double)(2*glist[empty]->inlinks) / (double)totallinks - */
/* 	    (double)((glist[empty]->totlinks + glist[empty]->inlinks) * */
/* 		     (glist[empty]->totlinks + glist[empty]->inlinks) ) / */
/* 	    (double)(totallinks * totallinks ); */

/* 	  dE += 2.0*(double)(glist[target]->inlinks+glist[empty]->inlinks+nlink) */
/* 	    / (double)totallinks - */
/* 	    (double)(glist[target]->totlinks + glist[target]->inlinks + */
/* 		     glist[empty]->totlinks + glist[empty]->inlinks ) * */
/* 	    (double)(glist[target]->totlinks + glist[target]->inlinks + */
/* 		     glist[empty]->totlinks + glist[empty]->inlinks ) / */
/* 	    (double)(totallinks*totallinks); */
	
/* 	  // Accept the change according to "inverse" Metroppolis. */
/* 	  // Inverse means that the algor is applied to the split */
/* 	  // and NOT to the merge! */
/* 	  if( (dE > 0.0) && ( prng_get_next(gen) > exp(-dE/T) ) ){ */
/* 	    MergeGroups(glist[target],glist[empty]); */
/* 	  } */
/* 	  else{ */
/* 	    energy -= dE; */
/* 	  } */

/* 	} */

/*       } */
/*     } */

/*     T = T * Ts; */
/*   } */

/* /\*   printf("energy = %lf\n",energy); *\/ */
  
/*   return CompressPart(part); */
/* } */


/* // merge = 0 => No group merging */
/* struct group *SACommunityIdentWeight(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int merge, struct prng *gen) */
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
/*   totallinks += CountLinksWeight(p); */
  
/*   for( i=1; i<nnod; i++ ) { */
/*     p = p->next; */

/*     nlist[i] = p; */
/*     trans[p->num] = i; */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*     AddNodeToGroup(glist[i],p); */
/*     totallinks += CountLinksWeight(p); */
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
/* 	target = floor(prng_get_next(gen) * nnod); */
/* 	g1 = nlist[target]->inGroup; */
	
/* 	if(glist[g1]->size < nnod){ */

/* 	  do{ */
/* 	    target = floor(prng_get_next(gen) * nnod); */
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
/* 	  if( (dE >= 0.0) || ( prng_get_next(gen) < exp(dE/T) ) ){ */
/* 	    MergeGroups(glist[g1],glist[g2]); */
/* 	    energy += dE; */
/* 	  } */
/* 	} */

/* 	// Split ///////////////////////////////////////////// */
/* 	target = floor(prng_get_next(gen) * nnod); // target node */
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
/* 	  if( (dE >= 0.0) && ( prng_get_next(gen) > exp(-dE/T) ) ){ */
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
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(prng_get_next(gen) * (double)nnod); */
/*       }while(newg == oldg); */

/*       // Calculate the change of energy */
/*       inold = CountLinksWeightInGroup(nlist[target],glist[oldg]); */
/*       innew = CountLinksWeightInGroup(nlist[target],glist[newg]); */
/*       nlink = CountLinksWeight(nlist[target]); */

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
/*       if( (dE >= 0.0) || ( prng_get_next(gen) < exp(dE/T) ) ){ */
/* 	MoveNode(nlist[target],glist[oldg],glist[newg]); */
/* 	energy += dE; */
/*       } */
/*     } */

/*     T = T * Ts; */
/*   } */

/* /\*   printf("energy = %lf\n",energy); *\/ */
  
/*   return CompressPart(part); */
/* } */


/* TranslationFromPart(struct group *part,int trans[]) */
/* { */
/*   struct group *g; */
/*   struct node_lis *node; */
/*   int count = 0; */
/*   int i; */

/*   for(i=0; i<maxim_int; i++){ */
/*     trans[i] = 0; */
/*   } */

/*   g = part; */
/*   while(g->next != NULL){ */
/*     g = g->next; */
/*     node = g->nodeList; */
/*     while(node->next != NULL){ */
/*       node = node->next; */

/*       trans[(node->ref)->num] = count; */
/*       count++; */
/*     } */
/*   } */
/* } */


/* double NodePR(struct node_gra *node, struct group *g) */
/* { */
/*   struct node_lis *nei = node->neig; */
/*   int ngroup = 0; */
/*   int ing[max_size]; */
/*   int nlink = 0; */
/*   double PR = 0.0; */
/*   int i; */

/*   for (i=0; i<max_size; i++){ */
/*     ing[i] = 0; */
/*   } */

/*   while(g->next != NULL){ */
/*     g = g->next; */
/*     if (g->size > 0) */
/*       ngroup++; */
/*   } */

/*   while(nei->next != NULL){ */
/*     nei = nei->next; */
    
/*     ing[nei->ref->inGroup] += 1; */
/*     nlink++; */
/*   } */

/*   if (nlink == 0) { */
/*     PR = 1.; */
/*   } */
/*   else { */
/*     for (i=0; i<max_size; i++){ */
/*       PR += (double)(ing[i] * ing[i]) / (double)(nlink * nlink); */
/*     } */
/*   } */

/* /\*   PR = (PR - 1.0/(double)ngroup) / (1.0 - 1.0/(double)ngroup); *\/ */

/*   return PR; */
/* } */


/* struct group *AllLevelsCommunities(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int merge, struct prng *gen) */
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
/* 	nettemp = BuildNetFromPart(g); */
	
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


/* InGroupBetZScore(struct group *part, double zsbet[]) */
/* { */
/*   struct node_gra *locnet,*p; */
/*   struct group *g; */
/*   double av,sigma; */
/*   int i; */
/*   double locbet[max_size]; */

/*   for (i=0; i<max_size; i++){ */
/*     zsbet[i] = 0.0; */
/*   } */

/*   g = part; */
/*   while(g->next != NULL){ */
/*     g = g->next; */
    
/*     // Build a network with the nodes in this community */
/*     locnet = BuildNetFromPart(g); */
    
/*     // Calculate the betweenness of the nodes in the community */
/*     CalculateBetweennessHistogram(locnet, locbet); */

/*     // Calculate average locbet */
/*     p = locnet; */
/*     av = 0.0; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       av += locbet[p->num]; */
/*     } */
/*     av /= (double)CountNodes(locnet); */

/*     // Calculate std deviation of locbet */
/*     p = locnet; */
/*     sigma = 0.0; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       sigma += (locbet[p->num] - av) * (locbet[p->num] - av); */
/*     } */
/*     sigma = sqrt(sigma / (double)CountNodes(locnet)); */
    
/*     // Calculate the zscore */
/*     p = locnet; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       zsbet[p->num] = (locbet[p->num] - av) / sigma; */
/*     } */

/*     RemoveGraph(locnet); */
/*   } */

/* } */


/* InGroupDegZScore(struct group *part, double z[]) */
/* { */
/*   struct node_gra *locnet=NULL, *p=NULL; */
/*   struct group *g=NULL; */
/*   double av=0.0, sigma=0.0; */
/*   int i; */

/*   g = part; */
/*   while(g->next != NULL){ */
/*     g = g->next; */
    
/*     locnet = BuildNetFromPart(g); */
    
/*     av = AverageDegree(locnet); */
/*     sigma = DegreeSTD(locnet); */
    
/*     p = locnet; */
/*     while((p = p->next) != NULL){ */
/*       if (sigma < 1.e-8) */
/* 	z[p->num] = 0.0; */
/*       else */
/* 	z[p->num] = ((double)CountLinks(p) - av) / sigma; */
/*     } */

/*     RemoveGraph(locnet); */
/*   } */
/* } */


/* InGroupDegFraction(struct group *part, double frac[]) */
/* { */
/*   struct node_gra *locnet,*p; */
/*   struct group *g; */
/*   int i; */

/*   g = part; */
/*   while(g->next != NULL){ */
/*     g = g->next; */
    
/*     locnet = BuildNetFromPart(g); */
    
/*     p = locnet; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
      
/*       frac[p->num] = (double)CountLinks(p) / (double)CountNodes(locnet); */
/*     } */

/*     RemoveGraph(locnet); */
/*   } */

/* } */

/* // Given a reference partition, calculates the fraction of */
/* // nodes that are correctly classified in an actual partition */
/* double CorrectlyClassified(struct group *refpart, */
/* 			   struct group *actpart) */
/* { */
/*   int i; */
/*   struct group *g1 = NULL; */
/*   struct group *g2 = NULL; */
/*   struct node_lis *p = NULL; */
/*   int group[maxim_int]; */
/*   int score[maxim_int]; */
/*   int map[maxim_int]; */
/*   int max_score[maxim_int]; */
/*   int correct = 0; */
/*   int nnod = 0; */
/*   double ratio; */

/*   for (i=0; i<maxim_int; i++){ */
/*     group[i] = -1; */
/*     map[i] = -1; */
/*     score[i] = 0; */
/*     max_score[i] = 0; */
/*   } */

/*   // Find the correct group of each node */
/*   g1 = refpart; */
/*   while(g1->next != NULL){ */
/*     g1 = g1->next; */

/*     p = g1->nodeList; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       group[p->node] = g1->label; */
/*     } */
/*   } */

/*   // Map the groups in the actual partition to the groups */
/*   // in the reference partition */
/*   g2 = actpart; */
/*   while(g2->next != NULL){ */
/*     g2 = g2->next; */
    
/*     // update the scores */
/*     p = g2->nodeList; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       score[group[p->node]] += 1; */
/*     } */

/*     // pick the group with a higher score */
/*     g1 = refpart; */
/*     while(g1->next != NULL){ */
/*       g1 = g1->next; */
/*       if( score[g1->label] > max_score[g2->label] ){ */
/* 	max_score[g2->label] = score[g1->label]; */
/* 	map[g2->label] = g1->label; */
/*       } */
/*       score[g1->label] = 0; */
/*     } */
/*   } */

/*   // If two groups are mapped into the same group of the  */
/*   // reference partition, get only the one that contains */
/*   // more nodes in that reference partition */
/*   g1 = actpart; */
/*   while(g1->next != NULL){ */
/*     g1 = g1->next; */
/*     g2 = g1; */
/*     while(g2->next != NULL){ */
/*       g2 = g2->next; */

/*       if (map[g1->label] == map[g2->label]){ */
/* 	if (max_score[g1->label] > max_score[g2->label]) */
/* 	  map[g2->label] = -1; */
/* 	else */
/* 	  map[g1->label] = -1; */
/*       } */
/*     } */
/*   } */

/*   // Count correctly classified nodes */
/*   g2 = actpart; */
/*   while(g2->next != NULL){ */
/*     g2 = g2->next; */
/*     nnod += g2->size; */

/*     if (map[g2->label] > -1){ */
/*       p = g2->nodeList; */
/*       while(p->next != NULL){ */
/* 	p = p->next; */
/* 	if(group[p->node] == map[g2->label]){ */
/* 	  correct++; */
/* 	} */
/* /\* 	else{ *\/ */
/* /\* 	  printf("%d missclassified\n",p->node+1); *\/ */
/* /\* 	} *\/ */
/*       } */
/*     } */
/* /\*     else{ *\/ */
/* /\*       p = g2->nodeList; *\/ */
/* /\*       while(p->next != NULL){ *\/ */
/* /\* 	p = p->next; *\/ */
/* /\* 	printf("%d missclassified\n",p->node+1); *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */

/*   } */

/*   ratio = (double)correct / (double)nnod; */

/*   return ratio; */
/* } */


/* struct group *KMeansRoleIdent(struct node_gra *net, int ngroup, int npass, char method, char distance, int bet_switch, int deg_switch, struct group *comm,int PR_switch, int z_switch, int locbet_switch, double *perf, double **cdata, double **csigma, int *dim) */
/* { */
/*   int i,j; */
/*   int nnod; */
/*   struct group *part = NULL; */
/*   struct group *glist[max_size]; */
/*   struct node_gra *p; */
/*   struct node_lis *ps; */
/*   int trans[maxim_int]; */
/*   int coun = 0; */
/*   double **nodprop; */
/*   int dsofar; */
/*   double bet[max_size]; */
/*   double locbet[max_size]; */
/*   double PR[max_size]; */
/*   double deg[max_size]; */
/*   double z[max_size]; */
/*   int degmin, degmax; */
/*   double betmin, betmax; */
/*   double locbetmin, locbetmax; */
/*   double PRmin, PRmax; */
/*   double zmin, zmax; */
/*   int **mask; */
/*   double weight[max_size]; */
/*   int result[max_size]; */
/*   double error; */
/*   int ifound; */

/*   *dim = 0; */

/*   nnod = CountNodes(net); */

/*   for(i=0; i<max_size; i++){ */
/*     bet[i] = -1.0; */
/*     locbet[i] = -1.0; */
/*     deg[i] = -1.0; */
/*     z[i] = -1.0; */
/*     PR[i] = -1.0; */
/*     weight[i] = 1.0; */
/*   } */

/*   // CHARACTERIZE THE DATA WITH THE SPECIFIED PARAMETERS */
/*   if (bet_switch != 0){ */
/*     *dim += 1; */
/*     CalculateBetweennessHistogram(net,bet); */
/*     betmin = 1e10; */
/*     betmax = -1e10; */
/*     p = net; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
      
/*       // CAUTION!!!! Log transformed betweenness */
/*       bet[p->num] = log(bet[p->num]); */

/*       if (bet[p->num] < betmin) betmin = bet[p->num]; */
/*       if (bet[p->num] > betmax) betmax = bet[p->num]; */
/*     } */
/*   } */

/*   if (deg_switch != 0){ */
/*     *dim += 1; */
/*     degmin = 1e6; */
/*     degmax = 0; */
/*     p = net; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       deg[p->num] = (double)CountLinks(p); */
/*       if (deg[p->num] < degmin) degmin = deg[p->num]; */
/*       if (deg[p->num] > degmax) degmax = deg[p->num]; */
/*     } */
/*   } */

/*   if (PR_switch != 0 && comm != NULL){ */
/*     *dim += 1; */
/*     p = net; */
/*     PRmin = 1.0; */
/*     PRmax = 0.0; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       PR[p->num] = NodePR(p,comm); */
/*       if (PR[p->num] < PRmin) PRmin = PR[p->num]; */
/*       if (PR[p->num] > PRmax) PRmax = PR[p->num]; */
/*     } */
/*   } */

/*   if (z_switch != 0 && comm != NULL){ */
/*     *dim += 1; */
/*     InGroupDegZScore(comm,z); */
/*     zmin = 1e10; */
/*     zmax = 0.0; */
/*     p = net; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       if (z[p->num] < zmin) zmin = z[p->num]; */
/*       if (z[p->num] > zmax) zmax = z[p->num]; */
/*     } */
/*   } */

/*   if (locbet_switch != 0 && comm != NULL){ */
/*     *dim += 1; */
/*     InGroupBetZScore(comm,locbet); */
/*     locbetmax = 1e100; */
/*     locbetmin = 0.0; */
/*     p = net; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       if (locbet[p->num] < locbetmin) locbetmin = locbet[p->num]; */
/*       if (locbet[p->num] > locbetmax) locbetmax = locbet[p->num]; */
/*     } */
/*   } */

/*   // Error message */
/*   if (*dim == 0){ */
/*     printf("Need at least one criterium to identify roles!!\n"); */
/*     return NULL; */
/*   } */

/*   // BUILD THE MATRIX WITH THE PROPERTIES OF THE NODES */
/*   // Initialize some matrices for the kcluster function */
/*   nodprop = allocate_d_mat(nnod,*dim); */
/*   mask = allocate_i_mat(nnod,*dim); */

/*   for(i=0; i<nnod; i++) */
/*     for(j=0; j<*dim; j++) */
/*       mask[i][j] = 1; */

/*   // Normalize the values and build the matrix */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     trans[p->num] = coun; */
/*     coun++; */

/*     dsofar = 0; */
/*     if (deg_switch != 0) */
/*       nodprop[trans[p->num]][dsofar++] = */
/* 	(double)(deg[p->num] - degmin) / (double)(degmax - degmin); */
/*     if (bet_switch != 0) */
/*       nodprop[trans[p->num]][dsofar++] = */
/* 	(bet[p->num] - betmin) / (betmax - betmin); */
/* /\*     if (PR_switch != 0 && comm != NULL) *\/ */
/* /\*       nodprop[trans[p->num]][dsofar++] = PR[p->num]; *\/ */
/*     if (PR_switch != 0 && comm != NULL) */
/*       nodprop[trans[p->num]][dsofar++] = */
/* 	(PR[p->num] - PRmin) / (PRmax - PRmin); */
/*     if (z_switch != 0 && comm != NULL) */
/*       nodprop[trans[p->num]][dsofar++] = (z[p->num] - zmin) / (zmax - zmin); */
/*     if (locbet_switch != 0 && comm != NULL) */
/*       nodprop[trans[p->num]][dsofar++] = (locbet[p->num] - locbetmin) / (locbetmax - locbetmin); */
/*   } */
/*   // The properties matrix has been built!! */

/*   // CREATE THE GROUPS AND ASSIGN THE NODES ACCORDING TO THE */
/*   // RESULT OF K-MEANS */
/*   kcluster(ngroup, nnod, *dim, nodprop, mask, weight, 0, npass, method, distance, result, cdata, &error, &ifound); */

/*   part = CreateHeaderGroup(); */
/*   ResetNetGroup(net); // All nodes reset to group -1 */

/*   // create the groups */
/*   glist[0] = CreateGroup(part,0); */
/*   for(i=1; i<ngroup; i++){ */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*   } */

/*   // place the nodes */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     AddNodeToGroup(glist[result[trans[p->num]]],p); */
/*   } */

/*   // Calculate the within group standard deviations */
/*   for(i=0; i<ngroup; i++){ */

/*     for(j=0; j<*dim; j++){ */
/*       csigma[i][j] = 0.0; */
/*     } */

/*     ps = glist[i]->nodeList; */
/*     while(ps->next != NULL){ */
/*       ps = ps->next; */

/*       for(j=0; j<*dim; j++){ */
/* 	csigma[i][j] += */
/* 	  (nodprop[trans[ps->node]][j] - cdata[i][j]) * */
/* 	  (nodprop[trans[ps->node]][j] - cdata[i][j]); */
/*       } */
/*     } */

/*     for(j=0; j<*dim; j++){ */
/*       csigma[i][j] = sqrt(csigma[i][j]/(double)glist[i]->size); */
/*     } */
/*   } */

/*   *perf = (double)ifound / (double)npass; */
  
/*   free_d_mat(nodprop,nnod); */
/*   free_i_mat(mask,nnod); */

/*   return CompressPart(part); */
/* } */


/* SOMRoleIdent(struct node_gra *net, int Nx, int Ny, int niter, char distance, int bet_switch, int deg_switch, struct group *comm,int PR_switch, int z_switch, int locbet_switch, double ***celldata, double **celldist, int clusterid[][2], int *dim, int inv_trans[]) */
/* { */
/*   int i, j, k; */
/*   int nnod; */
/*   struct group *part = NULL; */
/*   struct group *glist[max_size]; */
/*   struct node_gra *p; */
/*   struct node_lis *ps; */
/*   int coun = 0; */
/*   double **nodprop; */
/*   int dsofar; */
/*   double bet[max_size]; */
/*   double locbet[max_size]; */
/*   double PR[max_size]; */
/*   double deg[max_size]; */
/*   double z[max_size]; */
/*   int degmin, degmax; */
/*   double betmin, betmax; */
/*   double locbetmin, locbetmax; */
/*   double PRmin, PRmax; */
/*   double zmin, zmax; */
/*   int **mask; */
/*   double weight[max_size]; */
/*   double error; */
/*   int ifound; */
/*   int trans[maxim_int]; */
/*   int **normal; */
  
/*   *dim = 0; */

/*   nnod = CountNodes(net); */

/*   for(i=0; i<max_size; i++){ */
/*     bet[i] = -1.0; */
/*     locbet[i] = -1.0; */
/*     deg[i] = -1.0; */
/*     z[i] = -1.0; */
/*     PR[i] = -1.0; */
/*     weight[i] = 1.0; */
/*   } */

/*   // CHARACTERIZE THE DATA WITH THE SPECIFIED PARAMETERS */
/*   if (bet_switch != 0){ */
/*     *dim += 1; */
/*     CalculateBetweennessHistogram(net,bet); */
/*     betmin = 1e10; */
/*     betmax = -1e10; */
/*     p = net; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
      
/*       // CAUTION!!!! Log transformed betweenness */
/*       bet[p->num] = log(bet[p->num]); */

/*       if (bet[p->num] < betmin) betmin = bet[p->num]; */
/*       if (bet[p->num] > betmax) betmax = bet[p->num]; */
/*     } */
/*   } */

/*   if (deg_switch != 0){ */
/*     *dim += 1; */
/*     degmin = 1e6; */
/*     degmax = 0; */
/*     p = net; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       deg[p->num] = (double)CountLinks(p); */
/*       if (deg[p->num] < degmin) degmin = deg[p->num]; */
/*       if (deg[p->num] > degmax) degmax = deg[p->num]; */
/*     } */
/*   } */

/*   if (PR_switch != 0 && comm != NULL){ */
/*     *dim += 1; */
/*     p = net; */
/*     PRmin = 1.0; */
/*     PRmax = 0.0; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       PR[p->num] = NodePR(p,comm); */
/*       if (PR[p->num] < PRmin) PRmin = PR[p->num]; */
/*       if (PR[p->num] > PRmax) PRmax = PR[p->num]; */
/*     } */
/*   } */

/*   if (z_switch != 0 && comm != NULL){ */
/*     *dim += 1; */
/*     InGroupDegZScore(comm,z); */
/*     zmin = 1e10; */
/*     zmax = 0.0; */
/*     p = net; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       if (z[p->num] < zmin) zmin = z[p->num]; */
/*       if (z[p->num] > zmax) zmax = z[p->num]; */
/*     } */
/*   } */

/*   if (locbet_switch != 0 && comm != NULL){ */
/*     *dim += 1; */
/*     InGroupBetZScore(comm,locbet); */
/*     locbetmax = 0.0; */
/*     locbetmin = 1e100; */
/*     p = net; */
/*     while(p->next != NULL){ */
/*       p = p->next; */
/*       if (locbet[p->num] < locbetmin) locbetmin = locbet[p->num]; */
/*       if (locbet[p->num] > locbetmax) locbetmax = locbet[p->num]; */
/*     } */
/*   } */

/*   // Error message */
/*   if (*dim == 0){ */
/*     printf("Need at least one criterium to identify roles!!\n"); */
/* /\*     return NULL; *\/ */
/*   } */

/*   // BUILD THE MATRIX WITH THE PROPERTIES OF THE NODES */
/*   // Initialize some matrices for the somcluster function */
/*   nodprop = allocate_d_mat(nnod,*dim); */
/*   mask = allocate_i_mat(nnod,*dim); */
/*   normal = allocate_i_mat(Nx,Ny); */

/*   for(i=0; i<Nx; i++) */
/*     for(j=0; j<Ny; j++) */
/*       normal[i][j] = 0; */

/*   for(i=0; i<nnod; i++) */
/*     for(j=0; j<*dim; j++) */
/*       mask[i][j] = 1; */

/*   // Normalize the values and build the matrix */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     trans[p->num] = coun; */
/*     inv_trans[coun] = p->num; */
/*     coun++; */

/*     dsofar = 0; */
/*     if (deg_switch != 0) */
/*       nodprop[trans[p->num]][dsofar++] = */
/* 	(double)(deg[p->num] - degmin) / (double)(degmax - degmin); */
/*     if (bet_switch != 0) */
/*       nodprop[trans[p->num]][dsofar++] = */
/* 	(bet[p->num] - betmin) / (betmax - betmin); */
/* /\*     if (PR_switch != 0 && comm != NULL) *\/ */
/* /\*       nodprop[trans[p->num]][dsofar++] = PR[p->num]; *\/ */
/*     if (PR_switch != 0 && comm != NULL) */
/*       nodprop[trans[p->num]][dsofar++] = */
/* 	(PR[p->num] - PRmin) / (PRmax - PRmin); */
/*     if (z_switch != 0 && comm != NULL) */
/*       nodprop[trans[p->num]][dsofar++] = (z[p->num] - zmin) / (zmax - zmin); */
/*     if (locbet_switch != 0 && comm != NULL) */
/*       nodprop[trans[p->num]][dsofar++] = (locbet[p->num] - locbetmin) / (locbetmax - locbetmin); */

/*     printf("%d %d", p->num+1, inv_trans[coun-1]+1); */
/*     for (i=0; i<dsofar; i++) */
/*       printf(" %lf", nodprop[trans[p->num]][i]); */
/*     printf("\n"); */
/*   } */
/*   // The properties matrix has been built!! */

/*   // CREATE THE GROUPS AND ASSIGN THE NODES ACCORDING TO THE */
/*   // RESULT OF SOM */
/*   somcluster(nnod, *dim, nodprop, mask, weight, 0, Nx, Ny, 0.02, niter, distance, celldata, clusterid); */
  
/*   // Calculate the cell-distance matrix */
/*   for (i=0; i<Nx; i++){ */
/*     for (j=0; j<Ny; j++){ */
/*       celldist[i][j] = 0.0; */
/*       for (k=0; k<*dim; k++){ */
/* 	if (i > 0){ */
/* 	  celldist[i][j] += (celldata[i][j][k] - celldata[i-1][j][k]) */
/* 	    * (celldata[i][j][k] - celldata[i-1][j][k]); */
/* 	  normal[i][j] += 1; */
/* 	} */
/* 	if (i < Nx-1){ */
/* 	  celldist[i][j] += (celldata[i][j][k] - celldata[i+1][j][k]) */
/* 	    * (celldata[i][j][k] - celldata[i+1][j][k]); */
/* 	  normal[i][j] += 1; */
/* 	} */
/* 	if (j > 0){ */
/* 	  celldist[i][j] += (celldata[i][j][k] - celldata[i][j-1][k]) */
/* 	    * (celldata[i][j][k] - celldata[i][j-1][k]); */
/* 	  normal[i][j] += 1; */
/* 	} */
/* 	if (j < Ny-1){ */
/* 	  celldist[i][j] += (celldata[i][j][k] - celldata[i][j+1][k]) */
/* 	    * (celldata[i][j][k] - celldata[i][j+1][k]); */
/* 	  normal[i][j] += 1; */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   for (i=0; i<Nx; i++){ */
/*     for (j=0; j<Ny; j++){ */
/*       celldist[i][j] = sqrt(celldist[i][j] * (double)(*dim) / */
/* 			    (double)normal[i][j]);  */
/*     } */
/*   } */

/* /\*   part = CreateHeaderGroup(); *\/ */
/* /\*   ResetNetGroup(net); // All nodes reset to group -1 *\/ */

/* /\*   // create the groups *\/ */
/* /\*   glist[0] = CreateGroup(part,0); *\/ */
/* /\*   for(i=1; i<ngroup; i++){ *\/ */
/* /\*     glist[i] = CreateGroup(glist[i-1],i); *\/ */
/* /\*   } *\/ */

/* /\*   // place the nodes *\/ */
/* /\*   p = net; *\/ */
/* /\*   while(p->next != NULL){ *\/ */
/* /\*     p = p->next; *\/ */
/* /\*     AddNodeToGroup(glist[result[trans[p->num]]],p); *\/ */
/* /\*   } *\/ */

/* /\*   // Calculate the within group standard deviations *\/ */
/* /\*   for(i=0; i<ngroup; i++){ *\/ */

/* /\*     for(j=0; j<*dim; j++){ *\/ */
/* /\*       csigma[i][j] = 0.0; *\/ */
/* /\*     } *\/ */

/* /\*     ps = glist[i]->nodeList; *\/ */
/* /\*     while(ps->next != NULL){ *\/ */
/* /\*       ps = ps->next; *\/ */

/* /\*       for(j=0; j<*dim; j++){ *\/ */
/* /\* 	csigma[i][j] += *\/ */
/* /\* 	  (nodprop[trans[ps->node]][j] - cdata[i][j]) * *\/ */
/* /\* 	  (nodprop[trans[ps->node]][j] - cdata[i][j]); *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */

/* /\*     for(j=0; j<*dim; j++){ *\/ */
/* /\*       csigma[i][j] = sqrt(csigma[i][j]/(double)glist[i]->size); *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */

/* /\*   *perf = (double)ifound / (double)npass; *\/ */
  
/*   free_d_mat(nodprop,nnod); */
/*   free_i_mat(mask,nnod); */

/* /\*   return CompressPart(part); *\/ */
/* } */


/* struct group *CatalogRoleIdent(struct node_gra *net, struct group *comm) */
/* { */
/*   struct group *part = NULL; */
/*   struct node_gra *p; */
/*   double PR[maxim_int]; */
/*   double z[maxim_int]; */
/*   int ngroup = 7; */
/*   int dest_group; */
/*   int i; */
/*   struct group *glist[maxim_int]; */

/*   for(i=0; i<maxim_int; i++){ */
/*     z[i] = -1.0; */
/*     PR[i] = -1.0; */
/*     glist[i] = NULL; */
/*   } */

/*   // Characterize the nodes */
/*   InGroupDegZScore(comm, z); */
/*   p = net; */
/*   while((p = p->next) != NULL){ */
/*     PR[p->num] = NodePR(p,comm); */
/*     /\*     printf("%d %lf %lf\n", p->num+1, PR[p->num], z[p->num]); *\/ */
/*   } */

/*   // Create the roles and assign them according to the roles catalog */
/*   part = CreateHeaderGroup(); */
/*   ResetNetGroup(net); // All nodes reset to group -1 */

/*   // create the groups */
/*   glist[0] = CreateGroup(part,0); */
/*   for(i=1; i<ngroup; i++){ */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*   } */

/*   // place the nodes */
/*   p = net; */
/*   while((p = p->next) != NULL){ */
/*     if (z[p->num] <2.5){    // Node is not a hub */
/*       if (1.0 - PR[p->num] <= 0.05) */
/* 	dest_group = 0; */
/*       else if (1.0 - PR[p->num] < 0.6199) */
/* 	dest_group = 1; */
/*       else if (1.0 - PR[p->num] < 0.7999) */
/* 	dest_group = 2; */
/*       else */
/* 	dest_group = 3; */
/*     } */
/*     else{    // Node is a hub */
/*       if (1.0 - PR[p->num] < 0.2999) */
/* 	dest_group = 4; */
/*       else if (1.0 - PR[p->num] < 0.7499) */
/* 	dest_group = 5; */
/*       else */
/* 	dest_group = 6; */
/*     } */
    
/*     AddNodeToGroup(glist[dest_group], p); */
/*   } */

/*   return CompressPart(part); */
/* } */



/* struct group *PCatalogRoleIdent(struct node_gra *net, struct group *comm, int nbin) */
/* { */
/*   struct group *part = NULL; */
/*   struct node_gra *p; */
/*   double PR[maxim_int]; */
/*   double z[maxim_int]; */
/*   int ngroup = 100; */
/*   int dest_group; */
/*   int i; */
/*   struct group *glist[maxim_int]; */
/*   double bin = 1.0 / (double)nbin; */

/*   for(i=0; i<maxim_int; i++){ */
/*     z[i] = -1.0; */
/*     PR[i] = -1.0; */
/*     glist[i] = NULL; */
/*   } */

/*   // Characterize the nodes */
/*   InGroupDegZScore(comm,z); */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     PR[p->num] = NodePR(p,comm); */
/*   } */

/*   // Create the roles and assign them according to the roles catalog */
/*   part = CreateHeaderGroup(); */
/*   ResetNetGroup(net); // All nodes reset to group -1 */

/*   // create the groups */
/*   glist[0] = CreateGroup(part,0); */
/*   for(i=1; i<ngroup; i++){ */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*   } */

/*   // place the nodes */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */

/*     if (z[p->num] < 0){    // Node is an 'anti-hub' */
/*       dest_group = floor((1.0 - PR[p->num]) / bin); */
/*     } */
/*     else if (z[p->num] <2.5){    // Node is not a hub */
/*       dest_group = nbin + floor((1.0 - PR[p->num]) / bin); */
/*     } */
/*     else{    // Node is a hub */
/*       dest_group = 2 * nbin + floor((1.0 - PR[p->num]) / bin); */
/*     } */
    
/*     AddNodeToGroup(glist[dest_group],p); */
/*   } */

/*   return CompressPart(part); */
/* } */



/* RemoveInterGroupLinks(struct node_gra *net) */
/* { */
/*   struct node_gra *p; */
/*   struct node_lis *n; */
/*   int thisgroup; */

/*   p = net; */
  
/*   while(p->next != NULL){ */
/*     p = p->next; */
    
/*     thisgroup = p->inGroup; */

/*     n = p->neig; */
/*     while(n->next != NULL){ */
/*       n = n->next; */
      
/*       if(n->ref->inGroup != thisgroup) */
/* 	RemoveDirectedLink(p,n->ref); */
/*     } */
/*   } */
/* } */


/* struct group *ThermalNetworkSplitOdds(struct group *targ,double Ti, double Tf, struct prng *gen) */
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
/*   net = BuildNetFromPart(targ); */

/*   // Create the groups */
/*   split = CreateHeaderGroup(); */
/*   glist[0] = CreateGroup(split,0); */
/*   glist[1] = CreateGroup(split,1); */

/*   // Randomly assign the nodes to the groups */
/*   p = net; */
/*   while(p->next != NULL){ */
/*     p = p->next; */
/*     nlist[nnod] = p; */
/*     totallinks += CountLinks(p); */
/*     nnod++; */

/*     des = floor(prng_get_next(gen)*2.0); */
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
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       if(oldg == 0) */
/* 	newg = 1; */
/*       else */
/* 	newg = 0; */

/*       // Calculate the change of energy */
/*       inold = NLinksToGroup(nlist[target],glist[oldg]); */
/*       innew = NLinksToGroup(nlist[target],glist[newg]); */
/*       nlink = CountLinks(nlist[target]); */

/*       Anew = A - glist[newg]->size + innew + */
/* 	glist[oldg]->size - inold; */
/*       Bnew = B - innew + inold; */
/*       Cnew = C + glist[newg]->size - innew - */
/* 	glist[oldg]->size + inold; */
/*       Dnew = D + innew - inold;  */

/*       dE = Anew * Dnew / (Bnew * Cnew) - A * D / (B * C); */

/*       // Accept the change according to the Boltzman factor */
/*       if( (dE >= 0.0) || (prng_get_next(gen) < exp(dE/T)) ){ */
/* 	MoveNode(nlist[target],glist[oldg],glist[newg]); */
/*       } */
/*     } */

/*     T = T * Ts; */
/*   } */

/*   RemoveGraph(net); */

/*   return split; */
/* } */


/* // merge = 0 => No group merging */
/* struct group *SACommunityIdentOdds(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int merge, struct prng *gen) */
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
/*   totallinks += CountLinks(p); */
  
/*   for( i=1; i<nnod; i++ ) { */
/*     p = p->next; */

/*     nlist[i] = p; */
/*     trans[p->num] = i; */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*     AddNodeToGroup(glist[i],p); */
/*     totallinks += CountLinks(p); */
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
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(prng_get_next(gen) * (double)nnod); */
/*       }while(newg == oldg); */

/*       // Calculate the change of energy */
/*       inold = NLinksToGroup(nlist[target],glist[oldg]); */
/*       innew = NLinksToGroup(nlist[target],glist[newg]); */
/*       nlink = CountLinks(nlist[target]); */

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
/*       if( (dE >= 0.0) || ( prng_get_next(gen) < exp(dE/T) ) ){ */
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
/* 	target = floor(prng_get_next(gen) * (double)nnod); */
/* 	g1 = nlist[target]->inGroup; */
	
/* 	if(glist[g1]->size < nnod){ */

/* 	  do{ */
/* 	    target = floor(prng_get_next(gen) * (double)nnod); */
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
/* 	  if( (dE >= 0.0) || ( prng_get_next(gen) < exp(dE/T) ) ){ */
/* 	    MergeGroups(glist[g1],glist[g2]); */
/* 	    energy += dE; */
/* 	    A = Anew; */
/* 	    B = Bnew; */
/* 	    C = Cnew; */
/* 	    D = Dnew; */
/* 	  } */
/* 	} */

/* 	// Split ///////////////////////////////////////////// */
/* 	target = floor(prng_get_next(gen) * (double)nnod); // target node */
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
/* 	  if( (dE >= 0.0) && ( prng_get_next(gen) > exp(-dE/T) ) ){ */
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


/* struct group *SABlockModules(struct node_gra *net,double Ti,double Tf,double Ts,double fac,double alpha, struct prng *gen) */
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
/*   totallinks += CountLinks(p); */
  
/*   for( i=1; i<nnod; i++ ) { */
/*     p = p->next; */

/*     nlist[i] = p; */
/*     glist[i] = CreateGroup(glist[i-1],i); */
/*     AddNodeToGroup(glist[i],p); */
/*     totallinks += CountLinks(p); */
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
/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(prng_get_next(gen) * (double)nnod); */
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
/*       if( (dE >= 0.0) || ( prng_get_next(gen) < exp(dE/T) ) ){ */
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


/* struct group *SABlockModules3(struct node_gra *net,double Ti,double Tf,double Ts,double fac, int ngroup, struct prng *gen) */
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
/*     totallinks += CountLinks(p); */

/*     target = floor(prng_get_next(gen) * (double)ngroup); */
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

/*       target = floor(prng_get_next(gen) * (double)nnod); */
/*       oldg = nlist[target]->inGroup; */
/*       do{ */
/* 	newg = floor(prng_get_next(gen) * (double)ngroup); */
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
/*       if( (dE <= 0.0) || ( prng_get_next(gen) < exp(-dE/T) ) ){ */
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
/*       root_loc = BuildNetFromPart(g); */

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

/* PlotXfigModuleRole(struct node_gra *net, struct group *gpart, struct group *rpart, double mod_pos[][2], struct prng *gen) */
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
/* 	  0.5 * (double)(size * (prng_get_next(gen) - 0.5)); */
/* 	anode->coorY = ycent + */
/* 	  0.5 * (double)(size * (prng_get_next(gen) - 0.5)); */

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
	  
/* 	  target = floor(prng_get_next(gen) * (double)nspec); */
/* 	  agroup = specialg[target]; */
	
/* 	  do{ */
/* 	    dx = 0.01 * (double)(size * (prng_get_next(gen) - 0.5)); */
/* 	    dy = 0.01 * (double)(size * (prng_get_next(gen) - 0.5)); */
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
	  
/*       target = floor(prng_get_next(gen) * (double)nspec); */
/*       agroup = specialg[target]; */
      
/*       do{ */
/* 	dx = 0.01 * (double)(size * (prng_get_next(gen) - 0.5)); */
/* 	dy = 0.01 * (double)(size * (prng_get_next(gen) - 0.5)); */
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
/* struct node_gra *RandomizeNetworkRolePreserve(struct node_gra *net, struct group *part, double times, struct prng *gen) */
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
/*       niter =  ceil(times * (double)nlink + 0.5 + prng_get_next(gen)); */
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
/* 	    target1 = floor(prng_get_next(gen) * (double)nlink); */
/* 	    n1 = ori[target1]; */
/* 	    n2 = des[target1]; */
	   
/* /\* 	    printf("%d-%d\n", n1->num+1, n2->num+1); *\/ */
 
/* 	    target2 = floor(prng_get_next(gen) * (double)nlink); */
/* 	    if ( (g1 == g2) && (prng_get_next(gen) < 0.5) ) { */
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
/* 							  struct prng *gen) */
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
/*       niter =  ceil(times * (double)nlink + 0.5 + prng_get_next(gen)); */
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
/* 	      target1 = floor(prng_get_next(gen) * (double)nlink); */
/* 	      n1 = ori[target1]; */
/* 	      n2 = des[target1]; */
	   
/* /\* 	    printf("%d-%d\n", n1->num+1, n2->num+1); *\/ */
 
/* 	      target2 = floor(prng_get_next(gen) * (double)nlink); */
/* 	      if ( (g1 == g2) && (prng_get_next(gen) < 0.5) ) { */
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


/* double MutualInformation(struct group *part1, struct group *part2) */
/* { */
/*   struct group *g1 = NULL, *g2 = NULL; */
/*   struct node_lis *n1 = NULL, *n2 = NULL; */
/*   int S = 0, S1 = 0, S2 = 0, S12 = 0; */
/*   double H1 = 0.0, H2 = 0.0, H12 = 0.0; */
/*   double I12 = 0.0; */

/*   /\* */
/*     Count the number of nodes */
/*   *\/ */
/*   S = PartitionSize(part1); */
/*   if (S != PartitionSize(part2)) */
/*     printf("WARNING : partitions have different size!\n"); */

/*   /\* */
/*     Compute the H1 and H2 entropies */
/*   *\/ */
/*   g1 = part1; */
/*   while ((g1 = g1->next) != NULL) { */
/*     S1 = g1->size; */
/*     H1 += (double)S1 * log((double)S1 / (double)S); */
/*   } */
/*   g2 = part2; */
/*   while ((g2 = g2->next) != NULL) { */
/*     S2 = g2->size; */
/*     H2 += (double)S2 * log((double)S2 / (double)S); */
/*   } */

/*   /\* */
/*     Compute the join entropy H12 */
/*   *\/ */
/*   g1 = part1; */
/*   while ((g1 = g1->next) != NULL) { */
/*     S1 = g1->size; */

/*     g2 = part2; */
/*     while ((g2 = g2->next) != NULL) { */
/*       S2 = g2->size; */
      
/*       /\* */
/* 	Compute overlap */
/*       *\/ */
/*       S12 = 0; */
/*       n1 = g1->nodeList;  */
/*       while ((n1 = n1->next) != NULL) { */
/* 	n2 = g2->nodeList;  */
/* 	while ((n2 = n2->next) != NULL) { */
/* 	  if (n1->node == n2->node) { */
/* 	    S12++; */
/* 	    break; */
/* 	  } */
/* 	} */
/*       } // Overlap size S12 has been calculated */
/*       if (S12 > 0) */
/* 	H12 += (double)S12 * log((double)(S12 * S) / */
/* 				 (double)(S1 * S2)); */
/*     } */
/*   } // End of loop over groups */

/*   /\* */
/*     Compute mutual information */
/*   *\/ */
/*   I12 = -2.0 * H12 / (H1 + H2); */

/*   return I12; */
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
/* struct node_gra *SelectRandomNodeInGroup(struct group *mod, struct prng *gen) */
/* { */
/*   int i; */
/*   struct node_lis *node = mod->nodeList->next; */
/*   int target = floor(prng_get_next(gen) * mod->size); */

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


/* PrintPajekGraphPartition(struct node_gra *net) */
/* { */
/*   struct node_gra *p=net; */
/*   struct node_lis *li; */
/*   FILE *fit1,*fit2; */
/*   int trans[maxim_int]; */
/*   int transg[maxim_int]; */
/*   int i,co; */
/*   int nnod; */
/*   int gcount = 0; */

/*   for(i=0;i<maxim_int;i++){ */
/*     trans[i] = 0; */
/*     transg[i] = 0; */
/*   } */


/*   nnod = CountNodes(net); */

/*   fit1=fopen("graph.net","w"); */
/*   fit2=fopen("partition.clu","w"); */

/*   fprintf(fit1,"*Vertices %d\n",nnod); */
/*   fprintf(fit2,"*Vertices %d\n",nnod); */

/*   co=1; */
/*   p=net; */
/*   while(p->next!=NULL){ */
/*     p=p->next; */
/*     if(trans[p->num]==0){ */
/*       trans[p->num]=co; */
/*       co++; */
/*     } */
/*     fprintf(fit1,"%d *%d*           %lf %lf %lf\n",trans[p->num],p->num+1,p->coorX,p->coorY,p->coorZ); */
/* /\*     fprintf(fit1,"%d *%d*\n",trans[p->num],p->num+1); *\/ */

/*     if (transg[p->inGroup] == 0){ */
/*       gcount++; */
/*       transg[p->inGroup] = gcount; */
/*     } */
/* /\*     fprintf(fit2,"   %d\n",transg[p->inGroup]); *\/ */
/*     fprintf(fit2,"   %d\n", p->inGroup+1); */
/*   } */

/*   fprintf(fit1,"*Arcs\n"); */
/*   fprintf(fit1,"*Edges\n"); */
/*   p=net; */
/*   while(p->next!=NULL){ */
/*     p=p->next; */
/*     li=p->neig; */
/*     while(li->next!=NULL){ */
/*       li=li->next; */
/* /\*        if(p->num<li->node) *\/ */
/* 	fprintf(fit1,"%d   %d   1\n",trans[p->num],trans[li->node]); */
/*     } */
/*   } */

/*   fclose(fit1); */
/*   fclose(fit2); */
/* } */

/* PrintPajekGraphPartitionFName(struct node_gra *net, */
/* 			      char netf[], char parf[]) */
/* { */
/*   struct node_gra *p=net; */
/*   struct node_lis *li; */
/*   FILE *fit1,*fit2; */
/*   int trans[maxim_int]; */
/*   int i,co; */
/*   int nnod; */
/*   int gcount = 0; */

/*   for(i=0; i<maxim_int; i++){ */
/*     trans[i] = 0; */
/*   } */

/*   nnod = CountNodes(net); */

/*   fit1=fopen(netf,"w"); */
/*   fit2=fopen(parf,"w"); */

/*   fprintf(fit1,"*Vertices %d\n",nnod); */
/*   fprintf(fit2,"*Vertices %d\n",nnod); */

/*   co = 1; */
/*   p = net; */
/*   while ((p = p->next) != NULL) { */
/*     trans[p->num] = co++; */

/*     fprintf(fit1,"%d \"%d\"           %lf %lf %lf\n", */
/* 	    trans[p->num], p->num+1, p->coorX, p->coorY, p->coorZ); */
/* /\*     fprintf(fit1,"%d *%d*\n",trans[p->num],p->num+1); *\/ */

/*     fprintf(fit2,"   %d\n", p->inGroup+1); */
/*   } */

/*   fprintf(fit1,"*Arcs\n"); */
/*   fprintf(fit1,"*Edges\n"); */
/*   p=net; */
/*   while(p->next!=NULL){ */
/*     p=p->next; */
/*     li=p->neig; */
/*     while(li->next!=NULL){ */
/*       li=li->next; */
/* /\*        if(p->num<li->node) *\/ */
/*       fprintf(fit1,"%d   %d   1\n", */
/* 	      trans[p->num],trans[li->node]); */
/*     } */
/*   } */

/*   fclose(fit1); */
/*   fclose(fit2); */
/* } */

/* // Same as PrintPajekGraphPartition but the labels of the nodes */
/* // are taken from the "original" vector */
/* PrintPajekGraphPartitionTrans(struct node_gra *net,int label[]) */
/* { */
/*   struct node_gra *p=net; */
/*   struct node_lis *li; */
/*   FILE *fit1,*fit2; */
/*   int trans[max_size]; */
/*   int transg[max_size]; */
/*   int i,co; */
/*   int nnod; */
/*   int gcount = 0; */

/*   for(i=0;i<max_size;i++){ */
/*     trans[i] = 0; */
/*     transg[i] = 0; */
/*   } */

/*   nnod = CountNodes(net); */

/*   fit1=fopen("graph.net","w"); */
/*   fit2=fopen("partition.clu","w"); */

/*   fprintf(fit1,"*Vertices %d\n",nnod); */
/*   fprintf(fit2,"*Vertices %d\n",nnod); */

/*   co=1; */
/*   p=net; */
/*   while(p->next!=NULL){ */
/*     p=p->next; */
/*     if(trans[p->num]==0){ */
/*       trans[p->num]=co; */
/*       co++; */
/*     } */
/*     fprintf(fit1,"%d *%d*           %lf %lf %lf\n",trans[p->num],p->num+1,p->coorX,p->coorY,p->coorZ); */
/* /\*     fprintf(fit1,"%d *%d*\n",trans[p->num],label[p->num]+1); *\/ */

/*     if (transg[p->inGroup] == 0){ */
/*       gcount++; */
/*       transg[p->inGroup] = gcount; */
/*     } */
/*     fprintf(fit2,"   %d\n",transg[p->inGroup]); */
/*   } */

/*   fprintf(fit1,"*Arcs\n"); */
/*   fprintf(fit1,"*Edges\n"); */
/*   p=net; */
/*   while(p->next!=NULL){ */
/*     p=p->next; */
/*     li=p->neig; */
/*     while(li->next!=NULL){ */
/*       li=li->next; */
/* /\*        if(p->num<li->node) *\/ */
/* 	fprintf(fit1,"%d   %d   1\n",trans[p->num],trans[li->node]); */
/*     } */
/*   } */

/*   fclose(fit1); */
/*   fclose(fit2); */
/* } */


/* // Same as PrintPajekGraphPartition but groups smaller than */
/* // min are all marked as group 1 */
/* PrintPajekGraphPartitionMin(struct node_gra *net,struct group *part,int min) */
/* { */
/*   struct node_gra *p=net; */
/*   struct node_lis *li; */
/*   FILE *fit1,*fit2; */
/*   int trans[max_size]; */
/*   int transg[max_size]; */
/*   int i,co; */
/*   int nnod; */
/*   int gcount = 1; */
/*   struct group *glist[max_size]; */
/*   struct group *g = NULL; */

/*   for(i=0;i<max_size;i++){ */
/*     trans[i] = 0; */
/*     transg[i] = 0; */
/*     glist[i] = NULL; */
/*   } */

/*   g = part; */
/*   while(g->next != NULL){ */
/*     g = g->next; */
/*     glist[g->label] = g; */
/*   } */

/*   nnod = CountNodes(net); */

/*   fit1=fopen("graph.net","w"); */
/*   fit2=fopen("partition.clu","w"); */

/*   fprintf(fit1,"*Vertices %d\n",nnod); */
/*   fprintf(fit2,"*Vertices %d\n",nnod); */

/*   co=1; */
/*   p=net; */
/*   while(p->next!=NULL){ */
/*     p=p->next; */
/*     if(trans[p->num]==0){ */
/*       trans[p->num]=co; */
/*       co++; */
/*     } */
/*     fprintf(fit1,"%d *%d*           %lf %lf %lf\n",trans[p->num],p->num+1,p->coorX,p->coorY,p->coorZ); */
/* /\*     fprintf(fit1,"%d *%d*\n",trans[p->num],p->num+1); *\/ */

/*     if (transg[p->inGroup] == 0){ */
/*       if(glist[p->inGroup]->size >= min){ */
/* 	gcount++; */
/* 	transg[p->inGroup] = gcount; */
/*       } */
/*       else{ //groups smaller than min are pooled into partition 1 */
/* 	transg[p->inGroup] = 1; */
/*       } */
/*     } */
/*     fprintf(fit2,"   %d\n",transg[p->inGroup]); */
/*   } */

/*   fprintf(fit1,"*Arcs\n"); */
/*   fprintf(fit1,"*Edges\n"); */
/*   p=net; */
/*   while(p->next!=NULL){ */
/*     p=p->next; */
/*     li=p->neig; */
/*     while(li->next!=NULL){ */
/*       li=li->next; */
/* /\*        if(p->num<li->node) *\/ */
/* 	fprintf(fit1,"%d   %d   1\n",trans[p->num],trans[li->node]); */
/*     } */
/*   } */

/*   fclose(fit1); */
/*   fclose(fit2); */
/* } */
