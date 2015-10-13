#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_rng.h>

#include "tools.h"
#include "datastruct.h"
#include "graph.h"

// ---------------------------------------------------------------------
// A graph is a collection of nodes. This function creates an empty
// node to be placed at the top of the list.
// ---------------------------------------------------------------------
struct node_gra *
CreateHeaderGraph()
{
  struct node_gra *temp;

  temp = (struct node_gra *)calloc(1, sizeof(struct node_gra));
  temp->label = NULL;
  temp->num = -1;
  temp->coorX = -1.0;
  temp->coorY = -1.0;
  temp->coorZ = -1.0;
  temp->state = 0;
  temp->neig = NULL;
  temp->next = NULL;
  temp->inGroup = -1;
  temp->ivar1 = -1;
  temp->dvar1 = -1.;
  temp->trans = -1;
  temp->strength = 0;
  temp->degree = 0;

  return temp;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Node, link, and graph creation and memory allocation
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Create a node and add it at the end of the list that starts at p
// ---------------------------------------------------------------------
struct node_gra *
CreateNodeGraph(struct node_gra *p, char *label)
{

  while (p->next != NULL)
    p = p->next;
      
  p->next = (struct node_gra *) calloc(1, sizeof(struct node_gra));
  (p->next)->label = (char *) calloc(MAX_LABEL_LENGTH, sizeof(char));
  strcpy((p->next)->label, label);
  (p->next)->num = p->num + 1;
  (p->next)->state = 0;
  (p->next)->next = NULL;
  (p->next)->neig = (struct node_lis *) calloc(1, sizeof(struct node_lis));
  ((p->next)->neig)->node = -1;
  ((p->next)->neig)->next = NULL;
  (p->next)->ivar1 = 0;
  (p->next)->inGroup = -1;
  (p->next)->coorX = -1.0;
  (p->next)->coorY = -1.0;
  (p->next)->coorZ = -1.0;
  (p->next)->dvar1 = -1.0;
  (p->next)->trans = -1;
  (p->next)->strength = 0;
  (p->next)->degree = 0;

  return p->next;
}

// ---------------------------------------------------------------------
// The BFS list is a collection of node_bfs, used for breadth first
// searches. This function creates an empty node_bfs to be placed at
// the top of the list.
// ---------------------------------------------------------------------
struct node_bfs *
CreateHeaderList()
{
  struct node_bfs *temp;

  temp = (struct node_bfs *) calloc(1, sizeof(struct node_bfs));
  temp->d = -1;
  temp->ref = NULL;
  temp->next = NULL;
  temp->prev = NULL;
  temp->last = NULL;
  temp->pred = NULL;

  return temp;
}

// ---------------------------------------------------------------------
// Create an empty node for the tree data structure
// ---------------------------------------------------------------------
struct node_tree *
CreateNodeTree()
{
  struct node_tree *temp=NULL;
  temp = (struct node_tree *) calloc(1, sizeof(struct node_tree));
  temp->ref = NULL;
  temp->label = (char *) calloc(MAX_LABEL_LENGTH, sizeof(char));

  return temp;
}

// ---------------------------------------------------------------------
// Add node2 to the list of neighbors of node1: auto_link_sw=0
// prevents from self links being created; add_weight_sw!=0 enables
// the weight of the adjacency to increase if the adjacency already
// exists.
// ---------------------------------------------------------------------
int
AddAdjacency(struct node_gra *node1,
	     struct node_gra *node2,
	     int auto_link_sw,
	     int add_weight_sw,
	     double weight,
	     int status)
{
  struct node_lis *adja=NULL;

  // If node1==node2 and autolinks are not allowed, do nothing...
  if (node1 == node2 && auto_link_sw == 0) { 
    return 0;
  }
  // ...otherwise go ahead and try to create the link
  else {
	// Reset degree/strength of the nodes:
	node1->degree = 0;
	node2->degree = 0;
	node1->strength = 0;
	node2->strength = 0;
  	
    adja = node1->neig;
    while (adja->next != NULL) {
      adja = adja->next;
      if (adja->ref == node2) {
	// The link already exists
	if (add_weight_sw != 0) {
	  adja->weight += weight; // Update the weight of the link
	  return 1;
	}
	else {
	  return 0;            // Do nothing
	}
      }
    }

    // Create a new adjacency
    adja->next = (struct node_lis *)calloc(1, sizeof(struct node_lis));
    (adja->next)->node = node2->num;
    (adja->next)->nodeLabel = (char *)calloc(MAX_LABEL_LENGTH, sizeof(char));
    strcpy((adja->next)->nodeLabel, node2->label);
    (adja->next)->status = status;
    (adja->next)->next = NULL;
    (adja->next)->ref = node2;
    (adja->next)->btw = 0.0;
    (adja->next)->weight = weight;

    // Done
    return 1;
  }
}

// ---------------------------------------------------------------------
// Same as AddAdjacency, but we only pass the label of node2, and the
// ref pointer of the new adjacency is set to
// NULL. RewireAdjacencyByLabel needs to be run when AddAdjacencySoft
// is used
// ---------------------------------------------------------------------
int
AddAdjacencySoft(struct node_gra *node1,
		 char *node2Label,
		 int auto_link_sw,
		 int add_weight_sw,
		 double weight,
		 int status)
{
  struct node_lis *adja=NULL;

  // If node1==node2 and autolinks are not allowed, do nothing...
  if (strcmp(node1->label, node2Label) == 0 && auto_link_sw == 0) { 
    return 0;
  }
  // ...otherwise go ahead and try to create the adjacency
  else {
    adja = node1->neig;
    while (adja->next != NULL) {
      adja = adja->next;
      if (strcmp(adja->nodeLabel, node2Label) == 0) {
	// The link already exists
	if (add_weight_sw != 0) {
	  adja->weight += weight; // Update the weight of the link
	  return 1;
	}
	else {
	  return 0;            // Do nothing
	}
      }
    }
    
    // Create a new adjacency
    adja->next = (struct node_lis *) calloc(1, sizeof(struct node_lis));
    (adja->next)->nodeLabel = (char *) calloc(MAX_LABEL_LENGTH,
					      sizeof(char));
    strcpy((adja->next)->nodeLabel, node2Label);
    (adja->next)->status = status;
    (adja->next)->next = NULL;
    (adja->next)->ref = NULL;
    (adja->next)->btw = 0.0;
    (adja->next)->weight = weight;
    
    // Done
    return 1;
  }
}

// ---------------------------------------------------------------------
// Sets the ref pointers of all the adjacencies to the corresponding
// nodes based on their labels.
// ---------------------------------------------------------------------
void
RewireAdjacencyByLabel(struct node_gra *net)
{
  struct node_gra *p=NULL;
  struct node_lis *adja=NULL;
  void *nodeDict;

  // Map the nodes into a dictionary
  nodeDict = MakeLabelDict(net);

  // Point the adjacency pointers to the nodes
  p = net;
  while ((p = p->next) != NULL) {

	// Reset degree/strength of the nodes:
	p->degree = 0;
	p->strength = 0;
  
    adja = p->neig;
    while ((adja = adja->next) != NULL) {
      adja->ref = GetNodeDict(adja->nodeLabel, nodeDict);
      adja->node = adja->ref->num;
    }
  }

  // Done
  FreeLabelDict(nodeDict);
  return;
}

// ---------------------------------------------------------------------
// Sets the ref pointers of all the adjacencies to the corresponding
// nodes. You MUST use this after using AddAdjacencySoft!
// ---------------------------------------------------------------------
void
RewireAdjacencyByNum(struct node_gra *net)
{
  struct node_gra *p=NULL;
  struct node_lis *adja=NULL;
  struct node_gra **nlist;
  int nnod = CountNodes(net);

  // Map the nodes into an array for faster access
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  p = net;
  while ((p = p->next) != NULL) {
    nlist[p->num] = p;
  }

  // Point the adjacency pointers to the nodes and add the labels
  p = net;
  while ((p = p->next) != NULL) {
	// Reset degree/strength of the nodes:
	p->degree = 0;
	p->strength = 0;
  	
    adja = p->neig;
    while ((adja = adja->next) != NULL) {
      adja->ref = nlist[adja->node];
      strcpy(adja->nodeLabel, nlist[adja->node]->label);
    }
  }

  // Done
  free(nlist);
  return;
}

// ---------------------------------------------------------------------
// Takes the adjacency list of the origin node nori and copies it to
// the destination node ndes. ONLY soft links are created!!
// ---------------------------------------------------------------------
void
CopyAdjacencyList(struct node_gra *nori, struct node_gra *ndes)
{
  struct node_lis *pori=NULL;

  pori = nori->neig;
  while ((pori = pori->next) != NULL) {
    AddAdjacencySoft(ndes, pori->nodeLabel, 1, 0,
		     pori->weight, pori->status);
  }
}

// ---------------------------------------------------------------------
// Generates a copy of a given network.
// ---------------------------------------------------------------------
struct node_gra *
CopyNetwork(struct node_gra *p1)
{
  struct node_gra *root2=NULL, *p2=NULL;
  struct node_gra *last=NULL;

  // Copy the nodes and their adjacency lists
  root2 = CreateHeaderGraph();
  last = root2;
  while ((p1 = p1->next) !=  NULL) {
    p2 = CreateNodeGraph(last, p1->label);
    p2->coorX = p1->coorX;
    p2->coorY = p1->coorY;
    p2->coorZ = p1->coorZ;
    p2->state = p1->state;
    p2->ivar1 = p1->ivar1;
    p2->inGroup = p1->inGroup;
    p2->dvar1 = p1->dvar1;
    last = p2;
    CopyAdjacencyList(p1, p2);
  }

  /* Rewire the network to make all soft links hard and to set the
     labels */
  RewireAdjacencyByLabel(root2);

  /* Done */
  return root2;
}


// ---------------------------------------------------------------------
// Creates a binary tree for fast access to nodes by label
// ---------------------------------------------------------------------
void *
MakeLabelDict(struct node_gra *net)
{
  void *nodeDict = NULL;
  struct node_gra *p = net;
  struct node_tree *treeNode=NULL;

  while ((p = p->next) != NULL) {
    treeNode = CreateNodeTree();
    treeNode->label = strcpy(treeNode->label, p->label);
    treeNode = *(struct node_tree **)tsearch((void *)treeNode,
					     &nodeDict,
					     NodeTreeLabelCompare);
    treeNode->ref = p;
  }

  return nodeDict;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Node, link, and graph removal
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Frees the memory allocated to a node_tree.
// ---------------------------------------------------------------------
void
FreeNodeTree(struct node_tree *ntree)
{
  free(ntree->label);
  free(ntree);
  return;
}

// ---------------------------------------------------------------------
// Frees the memory allocated to a node_lis
// ---------------------------------------------------------------------
void
FreeNodeLis(struct node_lis *p)
{
  free(p->nodeLabel);
  free(p);
  return;
}

// ---------------------------------------------------------------------
// Recursively removes all the adjacencies in the adjacency list that
// hangs from p and frees the memory
// ---------------------------------------------------------------------
void
FreeAdjacencyList(struct node_lis *p)
{
  if (p->next != NULL) {
    FreeAdjacencyList(p->next);
  }
  FreeNodeLis(p);
}

// ---------------------------------------------------------------------
// Frees the memory allocated to a node_gra
// ---------------------------------------------------------------------
void
FreeNode(struct node_gra *node)
{
  if (node->neig != NULL) {
    FreeAdjacencyList(node->neig);
  }
  free(node->label);
  free(node);
  node = NULL;
  return;
}

// ---------------------------------------------------------------------
// Recursively removes all the nodes in a network and frees the memory
// ---------------------------------------------------------------------
void
RemoveGraph(struct node_gra *p)
{
  if (p->next != NULL) {
    RemoveGraph(p->next);
  }
  FreeNode(p);
}

// ---------------------------------------------------------------------
// Removes the link between nodes n1 and n2 (and frees the memory). If
// symmetric_sw != 0, RemoveLink will also remove the link from n2 to
// n1. CAUTION: RemoveLink will crash if there is no link between n1
// and n2 (or between n2 and n1 when symmetric_sw != 0).
// ---------------------------------------------------------------------
void
RemoveLink(struct node_gra *n1, struct node_gra *n2,
	   int symmetric_sw)
{
  struct node_lis *nn1;
  struct node_lis *nn2;
  struct node_lis *temp1;
  struct node_lis *temp2;

  // Reset degree/strength of the nodes:
  n1->degree = 0;
  n2->degree = 0;
  n1->strength = 0;
  n2->strength = 0;
  
  
  // Link n1-n2
  nn1 = n1->neig;
  while ((nn1->next)->ref != n2) {
    nn1 = nn1->next;
  }
  temp1 = nn1->next;
  nn1->next = temp1->next;
  FreeNodeLis(temp1);

  // Link n2-n1
  if (symmetric_sw != 0 && n1 != n2) {
    nn2 = n2->neig;
    while ((nn2->next)->ref != n1) {
      nn2 = nn2->next;
    }
    temp2 = nn2->next;
    nn2->next = temp2->next;
    FreeNodeLis(temp2);
  }
}

struct node_t
{
  /* Callers expect this to be the first element in the structure - do not
     move!  */
  const void *key;
  struct node_t *left;
  struct node_t *right;
  unsigned int red:1;
};

/**
@brief Free the memory allocated to a label dictionary.

It recursively walk the binary search tree and free its element.

@param dict A pointer to a dictionnary created with MakeLabelDict().
*/
void
FreeLabelDict(void *dict)
{
  struct node_t *focal = (struct node_t*) dict;
  if (focal != NULL){
	struct node_tree *leaf = (struct node_tree*) focal->key;
   	FreeLabelDict(focal->left); 
   	FreeLabelDict(focal->right);
	FreeNodeTree(leaf);
	free(focal);
  } 
}



// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Network input
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Compares two node_tree (required by tsearch)
// ---------------------------------------------------------------------
int
NodeTreeLabelCompare(const void *n1, const void *n2)
{
  return strcmp(((const struct node_tree *) n1)->label,
		((const struct node_tree *) n2)->label);
}

// ---------------------------------------------------------------------
// Builds a network from an imput file, which contains the list of
// links in the network (and, maybe, the weight of each link).
// ---------------------------------------------------------------------
struct node_gra *
FBuildNetwork(FILE *inFile,
	      int weight_sw,
	      int auto_link_sw,
	      int add_weight_sw,
	      int symmetric_sw)
{
  struct node_gra *root=NULL, *last_add=NULL;
  struct node_gra *n1=NULL, *n2=NULL;
  char label1[MAX_LABEL_LENGTH], label2[MAX_LABEL_LENGTH];
  void *node_dict=NULL;
  struct node_tree *n_tree=NULL, *ntree1=NULL, *ntree2=NULL;
  double weight;
  int noReadItems;
  char *line = NULL;
  size_t bufsiz = 0;
  ssize_t nbytes;

  // Create the header of the graph
  root = last_add = CreateHeaderGraph();
  
  // Go through the input file
  while ((nbytes = getline(&line, &bufsiz, inFile)) != -1){
    /* Read the labels (and weight, if necessary) */
    if (weight_sw == 0) {
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
   
    // Check if the nodes already exist, and create them otherwise
    n_tree = CreateNodeTree();
    strcpy(n_tree->label, label1);
    ntree1 = *(struct node_tree **)tsearch((void *)n_tree,
					   &node_dict,
					   NodeTreeLabelCompare);
    if (ntree1->ref == NULL) {
      ntree1->ref = last_add = CreateNodeGraph(last_add, label1);
    }
    else {
      FreeNodeTree(n_tree);
    }
    n1 = ntree1->ref;
    
    n_tree = CreateNodeTree();
    strcpy(n_tree->label, label2);
    ntree2 = *(struct node_tree **)tsearch((void *)n_tree,
					   &node_dict,
					   NodeTreeLabelCompare);
    if (ntree2->ref == NULL) {
      ntree2->ref = last_add = CreateNodeGraph(last_add, label2);
    }
    else {
      FreeNodeTree(n_tree);
    }
    n2 = ntree2->ref;
    
    // Add the link
    AddAdjacency(n1, n2, auto_link_sw, add_weight_sw, weight, 0);
    if (symmetric_sw != 0 && n1 != n2) {
      AddAdjacency(n2, n1, auto_link_sw, add_weight_sw, weight, 0);
    }
  }

  // Done
  free(line);
  FreeLabelDict(node_dict);
  return root;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Network printing and output
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Prints the degrees of all the nodes hanging from p
// ---------------------------------------------------------------------
void
FPrintDegrees(FILE *file, struct node_gra *p)
{
  while ((p = p->next) !=  NULL)
    fprintf(file, "%s %d\n", p->label, NodeDegree(p));
}

// ---------------------------------------------------------------------
// Prints a network as a list of adjacencies (pairs of nodes) with the
// corresponding weight if weight_sw is not 0. If symmetric_sw == 1,
// the link is only printed once (a-b), otherwise it will be printed
// twice if necessary (a-b, and b-a).
// ---------------------------------------------------------------------
void
FPrintNetAdjacencyList(FILE *outf,
		       struct node_gra *p,
		       int weight_sw,
		       int symmetric_sw)
{
  struct node_lis *n=NULL;
  int label_cmp;

  while((p = p->next) != NULL) {
    n = p->neig;
    while((n = n->next) != NULL) {
      label_cmp = strcmp(p->label, n->ref->label);
      if (weight_sw == 0) {
	if (symmetric_sw == 0 || label_cmp <= 0) {
	  fprintf(outf, "%s %s\n", p->label, n->ref->label);
	}
      }
      else {
	if (symmetric_sw == 0 || label_cmp <= 0) {
	  fprintf(outf, "%s %s %g\n", 
		  p->label, n->ref->label, n->weight);
	}
      }
    }
  }
}

/*
  ---------------------------------------------------------------------
  Prints a network in Pajek format. Node coordinates are printed if
  coor_sw is not 0. Weight is printed if weight_sw is not 0. If
  symmetric_sw == 1, the link is only printed once (a-b), otherwise
  it will be printed twice if necessary (a-b, and b-a).
  ---------------------------------------------------------------------
*/
void
FPrintPajekFile(char *fname,
		struct node_gra *root,
		int coor_sw,
		int weight_sw,
		int symmetric_sw)

{
  struct node_gra *p=root;
  struct node_lis *n=NULL;
  FILE *outF;

  outF = fopen(fname, "w");
  fprintf(outF, "*Vertices %d\n", CountNodes(root));

  /* Print the nodes */
  p = root;
  while((p = p->next) != NULL) {
    if (coor_sw == 1) {
      fprintf(outF, "%d \"%s\"           %lf %lf %lf\n",
	      p->num+1, p->label, p->coorX, p->coorY, p->coorZ);
    }
    else {
      fprintf(outF, "%d \"%s\"\n", p->num+1, p->label);
    }
  }
  
  /* Print the links title */
  if (symmetric_sw == 0) {
    fprintf(outF, "*Arcs\n");
  }
  else {
    fprintf(outF,"*Edges\n");
  }
  
  /* Print links */
  p = root;
  while((p = p->next) != NULL) {
    n = p->neig;
    while((n = n->next) != NULL) {
      if (weight_sw == 0) {
	if (symmetric_sw == 0 || p->num <= n->ref->num) {
	  fprintf(outF, "%d   %d\n", p->num+1, n->ref->num+1);
	}
      }
      else {
	if (symmetric_sw == 0 || p->num <= n->ref->num) {
	  fprintf(outF, "%d   %d   %g\n", 
		  p->num+1, n->ref->num+1, n->weight);
	}
      }
    }
  }

  /* Close the file and return */
  fclose(outF);
  return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Node, link, and graph operations
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Find and return a given node by number. The search is linear in the
// size of the network because it starts at the first node of the
// network and proceeds sequentially.
// ---------------------------------------------------------------------
struct node_gra *
GetNode(int num, struct node_gra *p)
{
  while((p->next)->num != num)
    p = p->next;
  return p->next;
}

// ---------------------------------------------------------------------
// Find and return a given node by label. Returns NULL if no node with
// that label is found. The search is logarithmic because a search
// tree is used.
// ---------------------------------------------------------------------
struct node_gra *
GetNodeDict(char *label, void *dict)
{
  struct node_tree *tempTreeNode = CreateNodeTree();
  void *treeNode = NULL;
  
  tempTreeNode->label = strcpy(tempTreeNode->label, label);
  treeNode = tfind((void *)tempTreeNode,
		   &dict,
		   NodeTreeLabelCompare);
  FreeNodeTree(tempTreeNode);
  
  if (treeNode != NULL)
    return (*(struct node_tree **)treeNode)->ref;
  else
    return NULL;
}

// ---------------------------------------------------------------------
// Find and return a given link
// ---------------------------------------------------------------------
struct node_lis *
GetLink(struct node_gra *n1, int n2)
{
  struct node_lis *pn = n1->neig;

  while((pn->next)->node !=  n2)
    pn = pn->next;
  return pn->next;
}

// ---------------------------------------------------------------------
// Returns 1 if the network hanging from p already contains a node
// with label "label"
// ---------------------------------------------------------------------
int
IsThereNode(char *label, struct node_gra *p)
{
  
  while ((p = p->next) != NULL)
    if (p->label == label)
      return 1;

  return 0;
}

/*
  -----------------------------------------------------------------------------
  Returns 1 if there is a link from node n1 to node n2
  -----------------------------------------------------------------------------
*/
int
IsThereLink(struct node_gra *n1, struct node_gra *n2)
{
  struct node_lis *n1n = n1->neig;

  while ((n1n = n1n->next) != NULL) {
    if (n1n->ref == n2) {
      return 1;
    }
  }

  return 0;
}

// ---------------------------------------------------------------------
// Returns 1 if there is a link between n1 and a node with number
// n2_num, and 0 otherwise
// ---------------------------------------------------------------------
int
IsThereLinkSoft(struct node_gra *n1, int n2_num)
{
  struct node_lis *pn;

  pn = n1->neig;
  while ((pn = pn->next) != NULL) {
    if(pn->node == n2_num)
      return 1;
  }

  return 0;
}

// ---------------------------------------------------------------------
// Removes isolated nodes from the network and returns the number of
// removed nodes. When done removing, nodes are renumbered so that
// they continue having consecutive integers. NEEDS TESTING.
// ---------------------------------------------------------------------
int
RemoveIsolatedNodes(struct node_gra *root)
{
  struct node_gra *p=NULL, *temp=NULL;
  int nrem=0;

  // Revome isolated nodes
  p = root;
  while (p->next !=  NULL) {
    if (NodeDegree(p->next) ==  0) {
      temp = p->next;
      p->next = (p->next)->next;
      free(temp);
      nrem++;
    }
    else{
      p = p->next;
    }
  }

  // Renumber the nodes and links in the network
  RenumberNodes(root);

  // Done
  return nrem;
}


// ---------------------------------------------------------------------
// Given a soft-linked network, remove all links that involve nodes
// that are not in the network.
// ---------------------------------------------------------------------
void
CleanAdjacencies(struct node_gra *net)
{
  struct node_lis *nei, *temp;
  struct node_gra *p = NULL;
  void *nodeDict;

  // Build a dictionary for fast access to nodes
  nodeDict = MakeLabelDict(net);

  // Clean the adjacencies
  p = net;
  while ((p = p->next) != NULL) {
	// Reset degree/strength of the nodes:
	p->degree = 0;
	p->strength = 0;
  
    nei = p->neig;
    while (nei->next != NULL) {
      if (GetNodeDict(nei->next->nodeLabel, nodeDict) == NULL) {
	temp = nei->next;
	if (temp->next == NULL)
	  nei->next = NULL;
	else
	  nei->next = temp->next;
	FreeNodeLis(temp);
      }
      else {
	nei = nei->next;
      }
    }
  }

  // Free memory
  FreeLabelDict(nodeDict);

  // Done
  return;
}

/*
  -----------------------------------------------------------------------------
  Remove nLinks random links in the network
  -----------------------------------------------------------------------------
*/
void
RemoveRandomLinks(struct node_gra *net,
		  int nLinks,
		  int symmetric_sw,
		  gsl_rng *gen)
{
  struct node_gra *p=net;
  struct node_lis *n=NULL;
  int *targetList, i, totalLinks, linkCounter=0;
  int n1;

  /* Create the list with the links to be removed */
  totalLinks = TotalNLinks(net, symmetric_sw);
  targetList = allocate_i_vec(totalLinks);
  for (i=0; i<totalLinks; i++) {
    targetList[i] = 0;
  }
  for (i=0; i<nLinks; i++) {
    do {
      n1 = gsl_rng_uniform_int(gen, totalLinks);
    } while (targetList[n1] == 1);
    targetList[n1] = 1;
  }
  
  /* Go over all links and remove those marked for removal. If links
     are undirected, we only try once per link (by imposing that the
     origin node has to be smaller than the destination node) */
  while ((p = p->next) != NULL) {
    n = p->neig;
    while (n->next != NULL) {
      if ((p->num <= ((n->next)->ref)->num) || (symmetric_sw == 0)) {
	if (targetList[linkCounter] == 1) {
	  RemoveLink(p, (n->next)->ref, symmetric_sw);
	}
	else {
	  n = n->next;
	}
	linkCounter++;
      }
      else {
	n = n->next;
      }
    }
  }

  /* Done */
  free_i_vec(targetList);
  return;
}

/*
  -----------------------------------------------------------------------------
  Add nLinks random links to the network
  -----------------------------------------------------------------------------
*/
void
AddRandomLinks(struct node_gra *net,
	       int nLinks,
	       int symmetric_sw,
	       gsl_rng *gen)
{
  int nnod, n;
  struct node_gra *p = net;
  struct node_gra **nlist;
  int n1, n2;

  /* Map the nodes to a list for faster access */
  nnod = CountNodes(net);
  nlist = (struct node_gra **) calloc(nnod, sizeof(struct node_gra *));
  p = net;
  while ((p = p->next) != NULL) {
    nlist[p->num] = p;
  }

  /* Add the links */
  for (n=0; n<nLinks; n++) {
    do {
      n1 = gsl_rng_uniform_int(gen, nnod);
      n2 = gsl_rng_uniform_int(gen, nnod);
    } while ((n1 == n2) || IsThereLink(nlist[n1], nlist[n2]) == 1);
    AddAdjacency(nlist[n1], nlist[n2], 0, 0, 1.0, 0);
    AddAdjacency(nlist[n2], nlist[n1], 0, 0, 1.0, 0);
  }

  /* Done */
  free(nlist);
  return;
}


/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  BFS list operations
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

// ---------------------------------------------------------------------
// Returns the node_bfs in list that corresponds to the node_gra node
// ---------------------------------------------------------------------
struct node_bfs *
GetBFS(struct node_gra *node, struct node_bfs *list)
{
  while (list != NULL) {
    if (list->ref == node) {
      return list;
    }
    else{
      list = list->next;
    }
  }

  return list;
}

// ---------------------------------------------------------------------
// ??????
// ---------------------------------------------------------------------
void
AddPredecessor(struct node_bfs *node, struct node_bfs *pred)
{
  struct pred *p = node->pred;
  
  while(p->next != NULL)
    p = p->next;

  p->next = (struct pred *)calloc(1,sizeof(struct pred));
  (p->next)->ref = pred;
  (p->next)->next = NULL;
}

// ---------------------------------------------------------------------
// Recursively remove all the predecessors in a predecessor list
// ---------------------------------------------------------------------
void
ClearPredecessors(struct pred *p)
{
  if (p->next != NULL)
    ClearPredecessors(p->next);
  free(p);
}

// ---------------------------------------------------------------------
// Count the number of predecessors of a given node in a BFS list
// ---------------------------------------------------------------------
int
CountPredecessors(struct node_bfs *node)
{
  struct pred *p = node->pred;
  int counter = 0;

  if((p->ref)->ref != NULL){  // Do not count the header as a
			      // predecessor
    while (p != NULL) {
      counter++;
      p = p->next;
    }
  }

  return counter;
}

// ---------------------------------------------------------------------
// Add a node_bfs to a list
// ---------------------------------------------------------------------
void
Enqueue(struct node_gra *node,
	struct node_bfs *predecessor,
	struct node_bfs *header,
	int *size,
	int dist)
{
  struct node_bfs *temp;
  
  if(node->state == 0){
    temp = (struct node_bfs *)calloc(1,sizeof(struct node_bfs));
    temp->d = dist;
    temp->next = NULL;
    temp->ref = node;
    temp->last = NULL;
    temp->pred = (struct pred *)calloc(1,sizeof(struct pred));
    (temp->pred)->ref = predecessor;
    (temp->pred)->next = NULL;
    *size += 1;
    node->state = 1;
    // Actualitzem l'apuntador a l'ultim element de la llista
    if(header->last == NULL){
      header->next = header->last = temp;
      temp->prev = header;
    }
    else{
      temp->prev = header->last;
      (header->last)->next = temp;
      header->last = temp;
    }
  }
  else{
    temp = GetBFS(node,predecessor);
    if((temp != NULL)&&(temp->d == dist))
      AddPredecessor(temp,predecessor);
  }
}

// ---------------------------------------------------------------------
// Add all the neighbors of a given node to a BFS list
// ---------------------------------------------------------------------
void
EnqueueAdjaList(struct node_bfs *lp,
		struct node_bfs *list,
		int *size,
		int d)
{
  struct node_lis *p=(lp->ref)->neig;

  while ((p = p->next) != NULL)
    Enqueue(p->ref,lp,list,size,d);
}

// ---------------------------------------------------------------------
// ????????
// ---------------------------------------------------------------------
struct node_gra *
DequeueOne(struct node_bfs *list,
	   struct node_bfs *one,
	   int *size)
{
  struct node_bfs *temp;
  struct node_bfs *p;
  struct node_gra *node;
  
  temp = one->next;

  if(list->last == temp){
    p = list;
    while(p->next != NULL){
      p->last = one;
      p = p->next;
    }
  }

  if(temp->next == NULL)
    one->next = NULL;
  else
    one->next = temp->next;

  *size -= 1;

  node = temp->ref;

  ClearPredecessors(temp->pred);
  free(temp);

  return node;
}

// ---------------------------------------------------------------------
// ???????
// ---------------------------------------------------------------------
struct node_gra *
Dequeue(struct node_bfs *list, int *size)
{
  struct node_bfs *temp;
  struct node_gra *node;

  temp = list->next;
  list->next = temp->next;
  if(list->next == NULL)
    list->last = NULL;
  else
    (list->next)->prev = list;

  *size -= 1;

  node = temp->ref;

  ClearPredecessors(temp->pred);
  free(temp);

  return node;
}

// ---------------------------------------------------------------------
// Update a BFS list with all the nodes at the following distance
// ---------------------------------------------------------------------
struct node_bfs *
RenewQueue(struct node_bfs *list,
	   struct node_bfs *lp,
	   int *size,
	   int d)
{
  while ((lp->next != NULL) && ((lp->next)->d == d)) {
    EnqueueAdjaList(lp->next, list, size, d+1);
    lp = lp->next;
  }
  return lp;
}

// ---------------------------------------------------------------------
// Dequeue all the nodes from a BFS list
// ---------------------------------------------------------------------
void
ClearList(struct node_bfs *list, int *size)
{
  struct node_gra *temp=NULL;

  while (list->next != NULL) {
    temp = Dequeue(list, size);
  }
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Network resetting
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Set the state of all nodes to 0
// ---------------------------------------------------------------------
void
ResetNodesState(struct node_gra *p)
{
  while ((p = p->next) != NULL)
    p->state = 0;
}

// ---------------------------------------------------------------------
// Make the numbers of the nodes consecutive
// ---------------------------------------------------------------------
void
RenumberNodes(struct node_gra *net)
{
  struct node_gra *p=NULL;
  struct node_lis *n=NULL;
  int count=0;

  // Renumber the nodes in the network
  p = net;
  while ((p = p->next) !=  NULL) {
    p->num = count++;
  }

  // Renumber the links (the adjacencies need to be pointed to the
  // reference nodes)
  p = net;
  while ((p = p->next) !=  NULL) {
    n = p->neig;
    while ((n = n->next) !=  NULL) {
      n->node = n->ref->num;
    }
  }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Network randomization
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

/*
  ---------------------------------------------------------------------
  Randomize the links of a network using the Markov chain switching
  algorithm
  ---------------------------------------------------------------------
*/
struct node_gra *
RandomizeSymmetricNetwork(struct node_gra *net,
			  double times,
			  gsl_rng *gen)
{
  int i;
  int nlink=0;
  int target1=0, target2=0;
  struct node_gra *p=NULL;
  struct node_lis *l=NULL;
  struct node_gra *n1, *n2, *n3, *n4;
  struct node_gra **ori, **des;
  int coun=0;
  int niter=0;

  /* Build the link lists (one for link origins and one for ends) */
  nlink = TotalNLinks(net, 1);
  niter =  ceil(times * (double)nlink + gsl_rng_uniform(gen));
  ori = (struct node_gra **)calloc(nlink, sizeof(struct node_gra *));
  des = (struct node_gra **)calloc(nlink, sizeof(struct node_gra *));
  p = net;
  while((p = p->next) !=  NULL){
    l = p->neig;
    while((l = l->next) !=  NULL){
      if (p->num > l->node) {
	ori[coun] = p;
	des[coun] = l->ref;
	coun++;
      }
    }
  }
  if(coun !=  nlink)
    fprintf(stderr, "Error in RandomizeNetwork: coun != nlink!!\n");

  /* Randomize the links */
  for (i=0; i<niter; i++) {
    /* select the 4 nodes */
    do {
      do {
	target1 = floor(gsl_rng_uniform(gen) * (double)nlink);
	n1 = ori[target1];
	n2 = des[target1];

	target2 = floor(gsl_rng_uniform(gen) * (double)nlink);
	if (gsl_rng_uniform(gen) < 0.5) {
	  n3 = des[target2];
	  n4 = ori[target2];
	}
	else {
	  n3 = ori[target2];
	  n4 = des[target2];
	}
      } while (n3 == n1 || n3 == n2 || n4 == n1 || n4 == n2);
    } while (IsThereLinkSoft(n1, n4->num) == 1 ||
	     IsThereLinkSoft(n2, n3->num) == 1);

    /* switch the link */
    RemoveLink(n1, n2, 1);
    RemoveLink(n3, n4, 1);
    AddAdjacency(n1, n4, 0, 0, 1., 1);
    AddAdjacency(n4, n1, 0, 0, 1., 1);
    AddAdjacency(n3, n2, 0, 0, 1., 1);
    AddAdjacency(n2, n3, 0, 0, 1., 1);

    ori[target1] = n1;
    des[target1] = n4;
    ori[target2] = n3;
    des[target2] = n2;
  }

  /* Free memory and return the network */
  free(ori);
  free(des);
  return net;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Network properties
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Counts the number of nodes in a network
// ---------------------------------------------------------------------
int
CountNodes(struct node_gra *p)
{
  int nodes = 0;
  while ((p = p->next) != NULL)
    nodes++;
  return nodes;
}



// ---------------------------------------------------------------------
// Get the degree of a node
// ---------------------------------------------------------------------

unsigned int NodeDegree(struct node_gra *node)
{
  if(!node->degree && node->neig->next != NULL){
    struct node_lis *p=node->neig;
    while ((p = p->next) != NULL)
      node->degree++;
  }
  return node->degree;
}


/*
  ---------------------------------------------------------------------
  Computes the average degree of the network, if symmetric_sw == 1, it
  will return the average in/out degree (NOT the average total
  degree!)
  ---------------------------------------------------------------------
*/
double
AverageDegree(struct node_gra *root, int symmetric_sw)
{
  if (symmetric_sw == 0)
    return (double)TotalNLinks(root, symmetric_sw) /
      (double)CountNodes(root);
  else
    return (double)(2 * TotalNLinks(root, symmetric_sw)) /
      (double)CountNodes(root);
}

/*
  ---------------------------------------------------------------------
  Computes the average of the squared degrees. If the network is
  directed, it will return the average of the square of the in/out
  degree (NOT the average total degree!)
  ---------------------------------------------------------------------
*/
double
AverageSquaredDegree(struct node_gra *root)
{
  struct node_gra *p=root;
  int k, nnod=0, k2sum=0;

  while ((p = p->next) != NULL) {
    nnod++;
    k = NodeDegree(p);
    k2sum += k * k;
  }

  return k2sum / nnod;

}

/*
   ---------------------------------------------------------------------
   Counts the number of links in the network. If symmetric_sw != 0,
   that is, if the network is symmetric, the number of links is
   divided by 2.
   ---------------------------------------------------------------------
*/
int
TotalNLinks(struct node_gra *p, int symmetric_sw)
{
  int total = 0;

  while((p = p->next) != NULL){
    total += NodeDegree(p);
  }
  
  if (symmetric_sw == 0)
    return total;
  else
    return total / 2;
}

/* 
   ---------------------------------------------------------------------
   Calculates the strength of a node, that is, the sum of the weights
   of all its links
   ---------------------------------------------------------------------
*/
double
NodeStrength(struct node_gra *node)
{
  if(!node->strength && node->neig->next != NULL){
    struct node_lis *p=node->neig;
    while ((p = p->next) != NULL)
      node->strength +=  p->weight;
  }
  return node->strength;
}

// ---------------------------------------------------------------------
// Print to a file the distribution of path lengths between all pairs
// of nodes. NEEDS TESTING.
// ---------------------------------------------------------------------
void
FPrintDistanceHistogram(FILE *file, struct node_gra *root)
{
  double histogram[10000];
  int a, d, nodes, *size, i, counter, size_ant;
  struct node_gra *p = root;
  struct node_bfs *list, *lp;
  double av_dis, norm;
  int **max_dis_pairs;

  max_dis_pairs = allocate_i_mat(CountNodes(root), 2);

  list = CreateHeaderList();

  for(i = 0;i<10000;i++)
    histogram[i] = 0.0;

  counter = 0;
  a = 0;
  size = &a;

  nodes = CountNodes(root);

  p = root;
  while(p->next != NULL){
    counter++;
    d = 0;
    printf("%d %d %d %d\n",counter,nodes,(p->next)->num,*size);
    ResetNodesState(root);
    p = p->next;

    Enqueue(p,list,list,size,0); //ATENCIO!! Posem el header com a
				 //predecessor del primer!!!!!
    lp = list;
    do{
      size_ant = *size;
      lp = RenewQueue(list,lp,size,d);
      d++;
      if(*size != size_ant){
	histogram[d] += (double)(*size-size_ant)/(double)nodes;
      }
    }while(*size != size_ant);

    ClearList(list, size);
  }

  // Output the histogram
  av_dis = norm = 0.0;
  for(i = 0;i<10000;i++){
    if(histogram[i] != 0.0){
      fprintf(file,"%d %lf\n",i,histogram[i]/(double)(nodes-1));
      av_dis += (double)(i)*histogram[i]/(double)(nodes-1);
      norm += (double)(histogram[i])/(double)(nodes-1);
    }
  }
  fprintf(file,"#Average distance = %lf\n",av_dis);
  fprintf(file,"#Normalization = %lf\n",norm);
  fprintf(file,"#z = %lf\n",histogram[1]);

  free_i_mat(max_dis_pairs, CountNodes(root));
  free(list);
}

// ---------------------------------------------------------------------
// Print to a file the distribution of path lengths between a node and
// all other nodes. NEEDS TESTING.
// ---------------------------------------------------------------------
void
FPrintDistanceHistogramFromNode(FILE *file,
				struct node_gra *root,
				int orinode)
{
  double histogram[10000];
  int a,d,nodes,*size,i,counter,size_ant;
  struct node_gra *p = root;
  struct node_bfs *list,*lp;
  double av_dis,norm;

  list = CreateHeaderList();

  for(i = 0;i<10000;i++)
    histogram[i] = 0.0;

  counter = 0;
  a = 0;
  size = &a;

  nodes = CountNodes(root);

  p = root;
  while(p->next != NULL){
    counter++;
    d = 0;
    p = p->next;

    if(p->num ==  orinode){
      ResetNodesState(root);
      
      Enqueue(p,list,list,size,0); //ATENCIO!! Posem el header com a
				   //predecessor del primer!!!!!
      lp = list;
      do{
	size_ant = *size;
	lp = RenewQueue(list,lp,size,d);
	d++;
	if(*size != size_ant){
	  histogram[d] += (double)(*size-size_ant);
	}
      }while(*size != size_ant);

      ClearList(list,size);
    }
  }

  av_dis = norm = 0.0;
  for(i = 0;i<10000;i++){
    if(histogram[i] != 0.0){
/*       fprintf(file,"%d %lf\n",i,histogram[i]); */
      av_dis += (double)(i)*histogram[i]/(double)(nodes-1);
      norm += (double)(histogram[i])/(double)(nodes-1);
    }
  }
  fprintf(file,"#Average distance = %lf\n",av_dis);
/*   fprintf(file,"#Normalization = %lf\n",norm); */
  fprintf(file,"#z = %lf\n",histogram[1]);

  free(list);
}

// ---------------------------------------------------------------------
// Calculates the average path length between nodes in a network
// ---------------------------------------------------------------------
double
AveragePathLength(struct node_gra *root)
{
  double histogram[10000];
  int a, d, nodes, *size, i, counter, size_ant;
  struct node_gra *p=root;
  struct node_bfs *list,*lp;
  double av_dis,norm;

  list = CreateHeaderList();

  for(i = 0;i<10000;i++)
    histogram[i] = 0.0;

  counter = 0;
  a = 0;
  size = &a;

  nodes = CountNodes(root);
  
  p = root;
  while(p->next != NULL){
    counter++;
    d = 0;
    ResetNodesState(root);
    p = p->next;

    Enqueue(p,list,list,size,0); //ATENCIO!! Posem el header com a
				 //predecessor del primer!!!!!
    lp = list;
    do{
      size_ant = *size;
      lp = RenewQueue(list,lp,size,d);
      d++;
      if(*size != size_ant){
	histogram[d] += (double)(*size-size_ant)/(double)nodes;
      }
    }while(*size != size_ant);

    if(*size != nodes){
      printf("Sorry, the network is not connected!\n");
      ClearList(list,size);
      free(list);
      return -1;
    }

    ClearList(list,size);
  }

  av_dis = norm = 0.0;
  for(i = 0;i<10000;i++){
    if(histogram[i] != 0.0){
      av_dis += (double)(i)*histogram[i]/(double)(nodes-1);
    }
  }

  free(list);

  return av_dis;
}

// ---------------------------------------------------------------------
// Calculates the average of the inverse of the path length between
// nodes in a network
// ---------------------------------------------------------------------
double
AverageInverseDistance(struct node_gra *root)
{
  double histogram[10000];
  int a,d,nodes,*size,i,counter,size_ant;
  struct node_gra *p = root;
  struct node_bfs *list,*lp;
  double av_dis,norm;

  list = CreateHeaderList();

  for(i = 0; i<10000; i++)
    histogram[i] = 0.0;

  counter = 0;
  a = 0;
  size = &a;

  nodes = CountNodes(root);
  
  p = root;
  while(p->next !=  NULL){
    counter++;
    d = 0;
    ResetNodesState(root);
    p = p->next;

    Enqueue(p,list,list,size,0); //ATENCIO!! Posem el header com a
				 //predecessor del primer!!!!!
    lp = list;
    do{
      size_ant = *size;
      lp = RenewQueue(list,lp,size,d);
      d++;
      if(*size !=  size_ant){
	histogram[d] +=  (double)(*size-size_ant)/(double)nodes;
      }
    }while(*size !=  size_ant);

    ClearList(list,size);
  }

  // Calculate the average inverse distance
  av_dis = norm = 0.0;
  for(i = 0; i<10000; i++){
    if(histogram[i] !=  0.0){
      av_dis +=  histogram[i]/(double)( i * (nodes-1) );
    }
  }

  free(list);

  return av_dis;
}


// ---------------------------------------------------------------------
// ??????????
// ---------------------------------------------------------------------
int
SumCommonLinks(struct node_gra *node, struct node_gra *root)
{
  int sum=0;
  struct node_lis *p;

  p = node->neig;
  while ((p = p->next) != NULL) {
    if ((p->ref)->state == 1) {
      sum++;
    }
  }

  return sum;
}

// ---------------------------------------------------------------------
// ??????????
// ---------------------------------------------------------------------
int CalculateLinksBetweenNeig(struct node_bfs *p,
			      struct node_gra *root)
{
  int sum = 0;

  while(p->next != NULL){
    p = p->next;
    sum += SumCommonLinks(p->ref,root);
  }

  return sum;
}

// ---------------------------------------------------------------------
// Calculates the clustering coefficent of a network as the average of
// the nodes individual clustering coefficient
// ---------------------------------------------------------------------
double
ClusteringCoefficient(struct node_gra *root)
{
  int nodes,*size,res_size;
  struct node_bfs *list;
  struct node_gra *p = root;
  double C_av;
  int k_v;
  int C_v;

  C_av = 0.0;
  res_size = 0;
  size = &res_size;
  list = CreateHeaderList();
  
  nodes = CountNodes(root);

  while(p->next != NULL){
    p = p->next;
    ResetNodesState(root);
    
    Enqueue(p,list,list,size,0); //Posem el header com a predecessor
				 //pero no afecta res!
    RenewQueue(list,list,size,0);
    Dequeue(list,size);
    k_v = *size;
    if(k_v>1){
      p->state = 0;
      C_v = CalculateLinksBetweenNeig(list,root);
      C_av += (double)C_v/(double)(k_v*(k_v-1));
    }
    else {
      nodes--;
    }
    ClearList(list,size);
  }

  C_av = C_av/(double)nodes;

  free(list);
  return C_av;
}

// ---------------------------------------------------------------------
// Calculates the clustering coefficent of a network as the total
// number of triangles in the network over the maximum number of
// possible triangles.
// ---------------------------------------------------------------------
double
ClusteringCoefficient2(struct node_gra *root)
{
  int nodes,*size,res_size;
  struct node_bfs *list;
  struct node_gra *p = root;
  int tri,ctr;  //#triangles and #connected triplets
  int k_v;
  double C;

  tri = ctr = 0;
  res_size = 0;
  size = &res_size;
  list = CreateHeaderList();
  
  nodes = CountNodes(root);

  while(p->next != NULL){
    p = p->next;
    ResetNodesState(root);
    
    Enqueue(p,list,list,size,0); //Posem el header com a predecessor
				 //pero no afecta res!
    RenewQueue(list,list,size,0);
    Dequeue(list,size);
    k_v = *size;
    if(k_v>1){
      p->state = 0;
      tri += CalculateLinksBetweenNeig(list,root);
      ctr += k_v*(k_v-1);
    }
    else
      nodes--;
    ClearList(list,size);
  }

  C = (double)tri/(double)ctr;

  free(list);
  return C;
}

// ---------------------------------------------------------------------
// Same as the clustering coefficient, but counting squares rather
// than triangles.
// ---------------------------------------------------------------------
double
SquareClustering(struct node_gra *root)
{
  int nodes,*size,res_size;
  struct node_bfs *list;
  struct node_bfs *lp;
  struct node_bfs *pl;
  struct node_gra *p = root;
  double C_av;
  int k1,k2;
  int C_v;
  int npre;
  int totnodes;
  int k1b;

  C_av = 0.0;
  res_size = 0;
  size = &res_size;
  list = CreateHeaderList();
  
  totnodes = nodes = CountNodes(root);

  while(p->next != NULL){
    p = p->next;
    ResetNodesState(root);
    
    Enqueue(p,list,list,size,0); //Posem el header com a predecessor
				 //pero no afecta res!
    lp = list;
    if(*size<totnodes)
      lp = RenewQueue(list,lp,size,0); // First neighbors
    else{ // If the size of the cluster is 1
      return -1.0;
    }
    k1b = *size-1;
    if(*size<totnodes){
      lp = RenewQueue(list,lp,size,1); // Second neighbors

      k1 = k2 = 0;
      pl = (list->next)->next; // Pointer to the first neighbor
      while(pl->d<2){
	pl = pl->next;
	k1++;
      }   // count first neighbors and place the pointer at the first
	  // second neighbor
      printf("%d %d = %d %d\n",p->num+1,k1,k1b,totnodes);
    
      C_v = 0;
      while(pl != NULL){
	k2++;
	npre = CountPredecessors(pl);
	C_v += npre*(npre-1);
	pl = pl->next;
      }

      if(k1>1){
	C_av += (double)(C_v)/(double)(k2*k1*(k1-1));
      }
      else
	nodes--;
    }
    else
      nodes--;

    ClearList(list,size);
  }

  C_av = C_av/(double)nodes;

  free(list);

  return C_av;
}

// ---------------------------------------------------------------------
// Compute the clustering coefficient of a single node
// ---------------------------------------------------------------------
double
OneNodeClusteringCoefficient(struct node_gra *node,
			     struct node_gra *root)
{
  int *size, res_size;
  struct node_bfs *list;
  int k_v;
  int C_v;
  double C;

  res_size = 0;
  size = &res_size;
  list = CreateHeaderList();
  
  ResetNodesState(root);
    
  Enqueue(node,list,list,size,0); // Posem el header com a predecessor
				  // pero no afecta res!
  RenewQueue(list,list,size,0);
  Dequeue(list,size);
  k_v = *size;
  if(k_v>1){
    node->state = 0;
    C_v = CalculateLinksBetweenNeig(list,root);
    C = (double)C_v/((double)k_v*(double)(k_v-1));
  }
  else
    C = -1.0000;

  ClearList(list,size);
  free(list);

  return C;
}

// ---------------------------------------------------------------------
// Compute the square clustering coefficient of a single node
// ---------------------------------------------------------------------
double
OneNodeSquareClustering(struct node_gra *node,
			struct node_gra *root)
{
  int *size, res_size;
  struct node_bfs *list;
  struct node_bfs *lp;
  struct node_bfs *pl;
  int k1,k2;
  int C_v;
  int npre;
  double T;
  int sizeant;

  res_size = 0;
  size = &res_size;
  list = CreateHeaderList();
  
  if(NodeDegree(node) == 0)
    return -1;

  ResetNodesState(root);
    
  Enqueue(node,list,list,size,0); //Posem el header com a predecessor
				  //pero no afecta res!
  lp = list;
  lp = RenewQueue(list,lp,size,0); // First neighbors

  sizeant = *size;

  lp = RenewQueue(list,lp,size,1); // Second neighbors

  if(*size !=  sizeant){

    k1 = k2 = 0;
    pl = (list->next)->next; // Pointer to the first neighbor
    while(pl->d<2){
      pl = pl->next;
      k1++;
    }   // count first neighbors and place the pointer at the first
	// second neighbor
    
    C_v = 0;
    while(pl != NULL){
      k2++;
      npre = CountPredecessors(pl);
      C_v += npre*(npre-1);
      pl = pl->next;
    }

    if(k1<= 1)
      return -2; // only 1 first neighbor!!
  }
  else
    return -3; // all nodes are first neighbors!!
  
  T = (double)C_v/((double)k2*(double)k1*(double)(k1-1));

  ClearList(list,size);
  free(list);

  return T;
}


/*
  ---------------------------------------------------------------------
  Calculates the betweenness of each link and stores in the btw field
  of the node_lis
  ---------------------------------------------------------------------
*/
void
CalculateLinkBetweenness(struct node_gra *root)
{
  int a,d,*size,i,counter,size_ant;
  struct node_gra *p = root;
  struct node_bfs *list,*lp;
  double *bet_loc;
  int n_pred;
  struct pred *ppre;
  int *ocupacio;
  struct node_lis *tar,*temp;
  int nnod=CountNodes(root);

  ocupacio = allocate_i_vec(nnod);
  bet_loc = allocate_d_vec(nnod);

  list = CreateHeaderList();

  for(i=0; i<nnod; i++){
    ocupacio[i] = 0;
  }
  counter = 0;
  a = 0;
  size = &a;

  // set the btw of all links to 0
  while(p->next != NULL){
    p = p->next;
    ocupacio[p->num] = 1;
    temp = p->neig;
    while(temp->next != NULL){
      temp = temp->next;
      temp->btw = 0.0;
    }
  }
  
  p = root;
  while(p->next != NULL){
    counter++;
    d = 0;
    ResetNodesState(root);
    p = p->next;

    Enqueue(p,list,list,size,0); //ATENCIO!! Posem el header com a
				 //predecessor del primer!!!!!
    lp = list;
    do{
      size_ant = *size;
      lp = RenewQueue(list,lp,size,d);
      d++;
    }while(*size != size_ant);

    //Update the betweenness
    for(i = 0;i<nnod;i++){
      bet_loc[i] = (double)ocupacio[i];
    }
    
    lp = list->last;
    while((n_pred = CountPredecessors(lp)) != 0){
      ppre = lp->pred;
      while(ppre != NULL){
	bet_loc[((ppre->ref)->ref)->num] +=
	  bet_loc[(lp->ref)->num]/(double)n_pred;
	
	tar = ((lp->ref)->neig)->next;
	while((tar->node) != ((ppre->ref)->ref)->num)
	  tar = tar->next;
	tar->btw += bet_loc[(lp->ref)->num]/(double)n_pred;

	tar = (((ppre->ref)->ref)->neig)->next;
	while((tar->node) != (lp->ref)->num)
	  tar = tar->next;
	tar->btw += bet_loc[(lp->ref)->num]/(double)n_pred;

	ppre = ppre->next;
      }
      lp = lp->prev;
    }
 
    ClearList(list,size);
  }
  
  free_i_vec(ocupacio);
  free_d_vec(bet_loc);
  free(list);
}

// ---------------------------------------------------------------------
// Calculates the betweenness of each link and writes in n1 and n2 the
// nodes corresponding to the largest betweenness in the network.
// ---------------------------------------------------------------------
void
CalculateBiggestLinkBetweenness(struct node_gra *root,int *n1,int *n2)
{
  int a,d,*size,i,counter,size_ant;
  struct node_gra *p = root;
  struct node_bfs *list,*lp;
  double *bet_loc;
  int n_pred;
  struct pred *ppre;
  int *ocupacio;
  struct node_lis *big = NULL;
  struct node_lis *tar,*temp;
  int nnod=CountNodes(root);

  ocupacio = allocate_i_vec(nnod);
  bet_loc = allocate_d_vec(nnod);
  list = CreateHeaderList();

  for(i = 0;i<nnod;i++){
    ocupacio[i] = 0;
  }
  counter = 0;
  a = 0;
  size = &a;

  while(p->next != NULL){
    p = p->next;
    ocupacio[p->num] = 1;
    temp = p->neig;
    while(temp->next != NULL){
      temp = temp->next;
      temp->btw = 0.0;
    }
  }
  
  big = ((root->next)->neig)->next;
  p = root;
  while(p->next != NULL){
    counter++;
    d = 0;
    ResetNodesState(root);
    p = p->next;

    Enqueue(p,list,list,size,0); //ATENCIO!! Posem el header com a
				 //predecessor del primer!!!!!
    lp = list;
    do{
      size_ant = *size;
      lp = RenewQueue(list,lp,size,d);
      d++;
    }while(*size !=  size_ant);

    //Actualize the betweenness
    for(i = 0;i<nnod;i++){
      bet_loc[i] = (double)ocupacio[i];
    }
    
    lp = list->last;
    while((n_pred = CountPredecessors(lp)) != 0){
      ppre = lp->pred;
      while(ppre != NULL){
	bet_loc[((ppre->ref)->ref)->num] +=
	  bet_loc[(lp->ref)->num]/(double)n_pred;
	
	tar = ((lp->ref)->neig)->next;
	while((tar->node) != ((ppre->ref)->ref)->num)
	  tar = tar->next;
	tar->btw += bet_loc[(lp->ref)->num]/(double)n_pred;
	if(big ==  NULL || tar->btw > big->btw){
	  big = tar;
	  *n2 = (lp->ref)->num;
	  *n1 = ((ppre->ref)->ref)->num;
	}

	tar = (((ppre->ref)->ref)->neig)->next;
	while((tar->node) != (lp->ref)->num)
	  tar = tar->next;
	tar->btw += bet_loc[(lp->ref)->num]/(double)n_pred;

	ppre = ppre->next;
      }
      lp = lp->prev;
    }
 
    ClearList(list,size);
  }

  free_i_vec(ocupacio);
  free_d_vec(bet_loc);
  free(list);
}

/*
  ----------------------------------------------------------------------------
  Calculates the betweenness of all nodes in a network. Results are
  stored in the dvar1 attribute of each node. The algorithm is
  directly taken from the pseudo-code in Algorithm 1 of the paper "A
  faster algorithm for betweenness centrality," by Ulrik Brandes.
  ----------------------------------------------------------------------------
*/
void
CalculateNodeBetweenness(struct node_gra *net)
{
  struct node_gra *nodep, *s, *v, *w;
  struct node_lis *wp;
  struct stack *S;
  struct queue *Q;
  struct queue **P;
  int *sigma, *d;
  double *delta;
  int nnod = CountNodes(net);
  int i;

  /* Reset dvar1 to 0 for all nodes */
  nodep = net;
  while ((nodep = nodep->next) != NULL)
    nodep->dvar1 = 0.0;

  /* Brandes algorithm */
  s = net;
  while ((s = s->next) != NULL) {
    S = stack_create();
    P = (struct queue **)calloc(nnod, sizeof(struct queue *));
    sigma = allocate_i_vec(nnod);
    d = allocate_i_vec(nnod);
    delta = allocate_d_vec(nnod);
    Q = queue_create();
    for (i=0; i<nnod; i++) {
      P[i] = queue_create();
      sigma[i] = 0;
      d[i] = -1;
      delta[i] = 0.0;
    }
    sigma[s->num] = 1;
    d[s->num] = 0;
    queue_enqueue(Q, (void *)s);
    while (!queue_empty(Q)) {
      v = (struct node_gra *)queue_dequeue(Q);
      stack_push(S, (void *)v);
      wp = v->neig;
      while ((wp = wp->next) != NULL) {
	w = wp->ref;
	/* w found for the first time? */
	if (d[w->num] < 0) {
	  queue_enqueue(Q, (void *)w);
	  d[w->num] = d[v->num] + 1;
	}
	/* Shortest path to w via v? */
	if (d[w->num] == d[v->num] + 1) {
	  sigma[w->num] = sigma[w->num] + sigma[v->num];
	  queue_enqueue(P[w->num], (void *)v);
	}
      }
    }

    /* S returns vertices in order of non-increasing distance from s */
    while (!stack_empty(S)) {
      w = (struct node_gra *)stack_pop(S);
      while (!queue_empty(P[w->num])) {
	v = (struct node_gra *)queue_dequeue(P[w->num]);
	delta[v->num] += (1. + delta[w->num]) * sigma[v->num] / sigma[w->num];
      }
      if (w->num != s->num)
	w->dvar1 += delta[w->num];
    }

    /* Free all memory */
    stack_free(S);
    for (i=0; i<nnod; i++) {
      queue_free(P[i]);
    }
    free(P);
    free_i_vec(sigma);
    free_i_vec(d);
    free_d_vec(delta);
    queue_free(Q);
  }
}

/*
  ---------------------------------------------------------------------
  Get some statistical properties (mean, std dev, min, and max) for
  the betweennesses of the nodes in the network.
  ---------------------------------------------------------------------
*/
void
NodeBetweennessStatistics(struct node_gra *net,
			  double *theMean,
			  double *theStddev,
			  double *theMin,
			  double *theMax)
{
  int nnod = CountNodes(net);
  double *betws = NULL;
  struct node_gra *p = net;
  
  /* Allocate memory */
  betws = allocate_d_vec(nnod);
  
  /* Get betweennesses */
  CalculateNodeBetweenness(net);
  while ((p = p->next) != NULL)
    betws[p->num] = p->dvar1;
  
  /* Get the statistical properties */
  (*theMean) = mean(betws, nnod);
  (*theStddev) = stddev(betws, nnod);
  (*theMin) = min(betws, nnod);
  (*theMax) = max(betws, nnod);

  /* Free memory and return */
  free_d_vec(betws);
  return;
}

// ---------------------------------------------------------------------
// Calculates the assortativity of a network
// ---------------------------------------------------------------------
double
Assortativity(struct node_gra *net)
{
  struct node_gra *p;
  struct node_lis *li;
  int *deg;
  int jk,jpk,j2pk2;
  int M;
  double r;

  // Initialize variables
  deg = allocate_i_vec(CountNodes(net));
  M = 0;
  jk = 0;
  jpk = 0;
  j2pk2 = 0;

  // Calculate nodes degrees
  p = net;
  while ((p = p->next) != NULL) {
    deg[p->num] = NodeDegree(p);
    M += deg[p->num];
  }

  // Calculate assortativity
  p = net;
  while ((p = p->next) != NULL) {
    li = p->neig;
    while ((li = li->next) != NULL) {
      jk += deg[p->num] * deg[(li->ref)->num];
      jpk += deg[p->num] + deg[(li->ref)->num];
      j2pk2 += deg[p->num] * deg[p->num] + 
	deg[(li->ref)->num] * deg[(li->ref)->num];
    }
  }

  r = ((double)jk / (double)M - ((double)jpk / (double)(2 * M)) *
       ((double)jpk / (double)(2 * M))) / 
    ((double)j2pk2 / (double)(2 * M) -
     ((double)jpk / (double)(2 * M)) *
     ((double)jpk / (double)(2 * M)));

  // Free memory and return
  free_i_vec(deg);
  return r;
}

// ---------------------------------------------------------------------
// Calculate the average degree of the neighbors of a node
// ---------------------------------------------------------------------
double
CalculateKnn(struct node_gra *node)
{
  struct node_lis *p=node->neig;
  int nneig=0, totdeg=0;

  while ((p = p->next) != NULL) {
    totdeg += NodeDegree(p->ref);
    nneig++;
  }
  
  return (double)totdeg / (double)nneig;
}

// ---------------------------------------------------------------------
// Returns 1 if the graph is connected and 0 if the graph has more
// than one component.
// ---------------------------------------------------------------------
int
IsGraphConnected(struct node_gra *p)
{
  struct node_bfs *list,*lp;
  int *size,r1;
  int size_ant;
  int d;
  int N = CountNodes(p);

  r1 = 0;
  size = &r1;

  ResetNodesState(p);
  p = p->next;
  list = CreateHeaderList();

  Enqueue(p,list,list,size,0);
  lp = list;
  d = 0;
  do{
    size_ant = *size;
    lp = RenewQueue(list,lp,size,d);
    d++;
  }while(*size != size_ant);

  ClearList(list,size);

  free(list);

  if(size_ant == N)
    return 1;
  else
    return 0;
}

// ---------------------------------------------------------------------
// ???????????
// ---------------------------------------------------------------------
int
AreConnectedList(struct node_gra *root,
		 struct node_gra *n1,
		 int cluslis[])
{
  struct node_bfs *list,*lp,*lp2;
  int *size,res1;
  int size_ant;
  int d;
  
  res1 = 0;
  size = &res1;

  list = CreateHeaderList();

  d = 0;
  ResetNodesState(root);
  Enqueue(n1,list,list,size,d); //ATENCIO!! Posem el header com a
				//predecessor del primer!!!!!
  lp = list;
  do{
    size_ant = *size;
    lp = RenewQueue(list,lp,size,d);
    d++;
    lp2 = lp;
    while(lp2->next != NULL){
      lp2 = lp2->next;
      if(cluslis[(lp2->ref)->num] == 1){
	ClearList(list,size);
	free(list);
	return 1;
      }
    }
  }while(*size != size_ant);
  
  ClearList(list,size);
  free(list);
  return 0;
}

// ---------------------------------------------------------------------
// Counts the number of strongly connected sets in the network
// ---------------------------------------------------------------------
int
CountStronglyConnectedSets(struct node_gra *root)
{
  struct node_gra *root_cop = NULL;
  struct node_gra *p;
  struct node_bfs *list,*lp,*lp2;
  int nnod = CountNodes(root);
  int anod = 0;
  int *size,res1;
  int *selected;
  int *this_net;
  int i;
  int size_ant;
  int d;
  int nclus = 0;

  selected = allocate_i_vec(nnod);
  this_net = allocate_i_vec(nnod);

  for(i = 0;i<nnod;i++)
    selected[i] = 0;
  
  root_cop = CopyNetwork(root);

  res1 = 0;
  size = &res1;

  list = CreateHeaderList();

  do{
    p = root->next;
    while(selected[p->num] != 0)
      p = p->next;
    selected[p->num] = 1;

    for(i = 0;i<nnod;i++)
      this_net[i] = 0;
    d = 0;

    ResetNodesState(root);
    Enqueue(p,list,list,size,d);
    anod++;
    lp = list;
    this_net[p->num] = 1;
    
    do{
      size_ant = *size;
      lp = RenewQueue(list,lp,size,d);
      lp2 = lp;
      while(lp2->next !=  NULL){
	if(AreConnectedList(root_cop,GetNode(((lp2->next)->ref)->num,root_cop),this_net) == 0){
	  DequeueOne(list,lp2,size);
	}
	else{
	  this_net[((lp2->next)->ref)->num] = 1;
	  selected[((lp2->next)->ref)->num] = 1;
	  anod++;
	  lp2 = lp2->next;
	}
      }
      d++;
    }while(*size !=  size_ant);

    ClearList(list,size);

    nclus++;

  }while(anod<nnod);

  free(list);
  RemoveGraph(root_cop);

  free_i_vec(selected);
  free_i_vec(this_net);
  return nclus;
}

/*
  ---------------------------------------------------------------------
  Creates a network that contains only the giant component of root
  ---------------------------------------------------------------------
*/
struct node_gra *
GetLargestStronglyConnectedSet(struct node_gra *root,
			       int thres)
{
  struct node_gra *root_cop=NULL;
  struct node_gra *root_loc=NULL;
  struct node_gra *p;
  struct node_gra *temp;
  struct node_bfs *list, *lp, *lp2;
  int nnod=CountNodes(root);
  int anod=0;
  int *size, res1;
  int *selected;
  int *this_net;
  int i;
  int size_ant;
  int d;
  int maxS = 0;
  struct node_gra *giant = NULL;

  selected = allocate_i_vec(nnod);
  this_net = allocate_i_vec(nnod);

  for(i = 0;i<nnod;i++)
    selected[i] = 0;
  
  root_cop = CopyNetwork(root);

  res1 = 0;
  size = &res1;

  list = CreateHeaderList();

  do{
    root_loc = CreateHeaderGraph();
    p = root->next;
    while(selected[p->num] != 0)
      p = p->next;
    selected[p->num] = 1;

    for(i = 0;i<nnod;i++)
      this_net[i] = 0;
    d = 0;

    ResetNodesState(root);
    Enqueue(p, list, list, size, d);
    anod++;
    lp = list;
    temp = CreateNodeGraph(root_loc, p->label);
    temp->num = p->num;
    CopyAdjacencyList(p, temp);
    this_net[p->num] = 1;
    
    do{
      size_ant = *size;
      lp = RenewQueue(list, lp, size, d);
      lp2 = lp;
      while(lp2->next != NULL){
	if(AreConnectedList(root_cop,
			    GetNode(((lp2->next)->ref)->num, root_cop),
			    this_net) == 0){
	  DequeueOne(list, lp2, size);
	}
	else{
	  temp = CreateNodeGraph(root_loc,((lp2->next)->ref)->label);
	  temp->num = lp2->next->ref->num;
	  CopyAdjacencyList(GetNode(((lp2->next)->ref)->num,root),temp);
	  this_net[temp->num] = 1;
	  selected[temp->num] = 1;
	  anod++;
	  lp2 = lp2->next;
	}
      }
      d++;
    }while(*size != size_ant);

    CleanAdjacencies(root_loc);
    RewireAdjacencyByLabel(root_loc);
    RenumberNodes(root_loc);

    if (CountNodes(root_loc) > maxS) {
      maxS = CountNodes(root_loc);
      if (giant !=  NULL) {
	RemoveGraph(giant);
      }
      giant = root_loc;
      if (maxS > nnod/2 || maxS > thres) {
	ClearList(list, size);
	RemoveGraph(root_cop); // WAS COMMENTED OUT?!?!?!?!
	return giant;
      }
    }
    else
      RemoveGraph(root_loc);

    ClearList(list,size);
    
  }while(anod<nnod);

  RemoveGraph(root_cop);

  free(list);
  free_i_vec(selected);
  free_i_vec(this_net);

  return giant;
}

struct node_gra *
GetLargestWeaklyConnectedSet(struct node_gra *root,int thres)
{
  struct node_gra *root_cop = NULL;
  struct node_gra *root_loc = NULL;
  struct node_gra *p, *p3;
  struct node_gra *temp;
  struct node_bfs *list,*lp,*lp2;
  int nnod=CountNodes(root);
  int anod=0;
  int *size, res1;
  int *selected;
  int *this_net;
  int i;
  int size_ant;
  int d;
  int maxS = 0;
  struct node_gra *giant = NULL;

  selected = allocate_i_vec(nnod);
  this_net = allocate_i_vec(nnod);

  for(i = 0;i<nnod;i++)
    selected[i] = 0;
  
  root_cop = CopyNetwork(root);

  res1 = 0;
  size = &res1;

  list = CreateHeaderList();
  
  // Reachable from p
  do{
    root_loc = CreateHeaderGraph();
    p = root->next;

    while(selected[p->num] != 0)
      p = p->next;
    selected[p->num] = 1;

    for(i = 0;i<nnod;i++)
      this_net[i] = 0;
    d = 0;

    ResetNodesState(root);
    Enqueue(p,list,list,size,d);
    anod++;
    lp = list;
    temp = CreateNodeGraph(root_loc, p->label);
    temp->num = p->num;
    CopyAdjacencyList(p,temp);
    this_net[p->num] = 1;
    
    do{
      size_ant = *size;
      lp = RenewQueue(list,lp,size,d);
      lp2 = lp;
      while(lp2->next != NULL){
	temp = CreateNodeGraph(root_loc,((lp2->next)->ref)->label);
	temp->num = lp2->next->ref->num;
	CopyAdjacencyList(GetNode(((lp2->next)->ref)->num,root),temp);
	this_net[temp->num] = 1;
	selected[temp->num] = 1;
	anod++;
	lp2 = lp2->next;
      }
      d++;
    }while(*size != size_ant);

    // Nodes that enable one to reach any node reachable from p
    p3 = p;

    while (p3->next != NULL){
      p3 = p3->next;

      if ((AreConnectedList(root_cop,GetNode(p3->num,root_cop),this_net) == 1)
	  && (selected[p3->num] == 0)) {

	temp = CreateNodeGraph(root_loc, p3->label);
	temp->num = p3->num;
	CopyAdjacencyList(p3, temp);
	this_net[temp->num] = 1;
	selected[temp->num] = 1;
	anod++;
	
	Enqueue(p3,list,list,size,d);

	do{
	  size_ant = *size;
	  lp = RenewQueue(list,lp,size,d);
	  lp2 = lp;
	  while(lp2->next != NULL){
	    temp = CreateNodeGraph(root_loc,((lp2->next)->ref)->label);
	    temp->num = lp2->next->ref->num;
	    CopyAdjacencyList(GetNode(((lp2->next)->ref)->num,root),temp);
	    this_net[temp->num] = 1;
	    selected[temp->num] = 1;
	    anod++;
	    lp2 = lp2->next;
	  }
	  d++;
	}while(*size != size_ant);

      }
    }

    // Rewire the cluster
    CleanAdjacencies(root_loc);
    RewireAdjacencyByLabel(root_loc);
    RenumberNodes(root_loc);

    // Check if this is the largest component
    if(CountNodes(root_loc)>maxS){
      maxS = CountNodes(root_loc);
      if (giant !=  NULL){
	RemoveGraph(giant);
      }
      giant = root_loc;
      if(maxS>nnod/2 || maxS>thres ){
	ClearList(list,size);
/* 	RemoveGraph(root_cop); */
	return giant;
      }
    }
    else
      RemoveGraph(root_loc);

    ClearList(list,size);

  }while(anod<nnod);

  RemoveGraph(root_cop);

  free(list);
  free_i_vec(selected);
  free_i_vec(this_net);
  return giant;
}

//Function to remove those links that go to other components
////It is necessary for the next function
void ClearAdjacencies(struct node_gra *p,int *net)
{
  struct node_lis *nei,*temp;

  while(p->next!=NULL){
    p=p->next;
    nei=p->neig;
    while(nei->next!=NULL){
      if(net[(nei->next)->node]==0){
        temp=nei->next;
        if(temp->next==NULL)  nei->next=NULL;
        else  nei->next=temp->next;
        free(temp);
      }else{
        nei=nei->next;
      }
    }
  }
}


// Returns the number of connected components found. 
// net_list must have been initialized, eg. struct node_gra *llista[1000];
int GetAllConnectedSets(struct node_gra *network,struct node_gra **net_list)
{
  struct node_gra *network_cop=NULL;
  struct node_gra *network_loc=NULL;
  struct node_gra *p;
  struct node_gra *temp, *temp2;  
  struct node_bfs *list,*lp,*lp2;
  int numNodes=0;
  int anod=0;
  int *size,res1;
  int *selected;
  int *thisNet;
  int i;
  int sizeAnt;
  int d;

  int numComponents=0;
  int MAX_SIZE_NET=10001;

  selected=(int *)malloc(sizeof(int)*MAX_SIZE_NET);
  thisNet=(int *)malloc(sizeof(int)*MAX_SIZE_NET);

  for(i=0;i<MAX_SIZE_NET;i++)
    selected[i]=0;

  network_cop=CopyNetwork(network);

  res1=0;
  size=&res1;

  numNodes=CountNodes(network);
  //printf("Network has %i nodes\nDetecting connected components\n",numNodes);
  list=CreateHeaderList();



  do{
    network_loc=CreateHeaderGraph();
    p=network->next;
    while(selected[p->num]!=0)  //Busquem 1 node no tractat encara
      p=p->next;
    selected[p->num]=1;


    for(i=0;i<MAX_SIZE_NET;i++)
      thisNet[i]=0;
    d=0;

    ResetNodesState(network);
    Enqueue(p,list,list,size,d);
    anod++;
    //printf("Starting network with node %d    Total: %d/%d\n",p->num,anod,numNodes);

    lp=list;
    temp=CreateNodeGraph(network_loc,p->label);
    temp->num=p->num;
    temp->state=p->state;
    CopyAdjacencyList(p,temp);
    thisNet[p->num]=1;

    do{

      sizeAnt=*size;
      lp=RenewQueue(list,lp,size,d);

      lp2=lp;
      while(lp2->next!=NULL){
        if(AreConnectedList(network_cop,GetNode(((lp2->next)->ref)->num,network_cop),thisNet)==0){
          DequeueOne(list,lp2,size);
        }else{
          temp=CreateNodeGraph(network_loc,((lp2->next)->ref)->label);
          temp2=GetNode(((lp2->next)->ref)->num,network);
          temp->num=temp2->num;
          temp->state=temp2->state;
          CopyAdjacencyList(temp2,temp);
          thisNet[temp->num]=1;
          selected[temp->num]=1;
          anod++;
          //printf("   adding %d    Total: %d/%d\n",temp->num,anod,numNodes);
          lp2=lp2->next;
        }
      }
      d++;
    }while(*size!=sizeAnt);

    ClearAdjacencies(network_loc,thisNet);
    RewireAdjacencyByLabel(network_loc);

    net_list[numComponents]=network_loc;

    //printf("Component %3i: Size %i\n",numComponents+1,CountNodes(net_list[numComponents]));
    numComponents++;
    ClearList(list,size);
  }while(anod<numNodes);

  free(list);
  RemoveGraph(network_cop);

  return numComponents;

}



/*
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Node and network comparison
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
*/

/*
  ---------------------------------------------------------------------
  Calculates the number of common neighbors between two nodes.
  ---------------------------------------------------------------------
*/
int
CommonNeighbors(struct node_gra *n1, struct node_gra *n2)
{
  int ncom=0;
  struct node_lis *p;

  p = n1->neig;
  while ((p = p->next) != NULL)
    ncom += IsThereLink(p->ref, n2);

  return ncom;
}

/*
  ---------------------------------------------------------------------
  Calculates Jaccard index between two nodes.
  ---------------------------------------------------------------------
*/
double
JaccardIndex(struct node_gra *n1, struct node_gra *n2)
{
  int ncom=0;
  struct node_lis *p;
  int k1, k2;

  /* Degrees of both nodes */
  k1 = NodeDegree(n1);
  k2 = NodeDegree(n2);

  /* Count the number of common neighbors */
  p = n1->neig;
  while ((p = p->next) != NULL) {
    if (p->ref == n2) {
      k1--;
      k2--;
    }
    else {
      ncom += IsThereLink(p->ref, n2);
    }
  }

  /* Evaluate Jaccard index */
  if (ncom == 0)
    return 0.0;
  else
    return (double)ncom / (double)(k1 + k2 - ncom);
}

/*
  ---------------------------------------------------------------------
  Calculates the topological overlap, as defined by Ravasz et al.,
  between two nodes.
  ---------------------------------------------------------------------
*/
double
TopologicalOverlap(struct node_gra *n1, struct node_gra *n2)
{
  double overlap;
  int ncom=0;
  int minim;
  struct node_gra *node;
  struct node_lis *p;
  int k1, k2;

  /* Degrees of both nodes */
  k1 = NodeDegree(n1);
  k2 = NodeDegree(n2);

  /* Determine the node with less neighbors */
  if (k1 < k2) {
    minim = k1;
    p = n1->neig;
    node = n2;
  }
  else{
    minim = k2;
    p = n2->neig;
    node = n1;
  }

  /* Count the number of common neighbors */
  while (p->next != NULL) {
    p = p->next;
    if (p->ref == node)
      ncom++;
    else
      ncom += IsThereLinkSoft( node, p->node);
  }

  /* Evaluate overlap */
  overlap = (double)ncom / (double)minim;

  return overlap;
}

/*
  -----------------------------------------------------------------------------
  Given two networks whose nodes are comparable (that is, their labels
  are the same in both networks), the subroutine returns a number of
  conditional probabilities.
  -----------------------------------------------------------------------------
*/
void
CompareTwoNetworks(struct node_gra *netA, struct node_gra *netB,
		   int *nA_nod, int *nB_nod, int *ncom_nod, 
		   double *p_nAB, double *p_nBA,
		   int *nA_lin, int *nB_lin, int *ncom_lin, 
		   double *p_lAB, double *p_lBA)
{
  struct node_gra *pA, *pB, *n1, *n2;
  struct node_lis *lA, *lB;
  int ntot_lin = 0;
  void *node_dictA=NULL, *node_dictB=NULL;
  struct node_tree *n_tree=NULL, *ntreeA=NULL, *ntreeB=NULL;

  /* Initialize variables */
  *ncom_nod = 0;
  *nA_nod = 0;
  *nB_nod = 0;
  *ncom_lin = 0;
  *nA_lin = 0;
  *nB_lin = 0;

  /*
    Calculate node overlap
  */
  /* Go through the first network and build its dictionary */
  pA = netA;
  while ((pA = pA->next) != NULL) {
    n_tree = CreateNodeTree();
    strcpy(n_tree->label, pA->label);
    ntreeA = *(struct node_tree **)tsearch((void *)n_tree,
					   &node_dictA,
					   NodeTreeLabelCompare);
    ntreeA->ref = pA;
    (*nA_nod) += 1;
  }
  
  /* Go through the second network and build its
     dictionary. Additionally, keep counting the number of common
     nodes (that is, how many nodes in B are also in A). */
  pB = netB;
  while ((pB = pB->next) != NULL) {
    n_tree = CreateNodeTree();
    strcpy(n_tree->label, pB->label);
    ntreeB = *(struct node_tree **)tsearch((void *)n_tree,
					   &node_dictB,
					   NodeTreeLabelCompare);
    ntreeB->ref = pB;
    (*nB_nod) += 1;

    /* Is this node in A? */
    if (GetNodeDict(pB->label, node_dictA) != NULL) {
      (*ncom_nod) += 1;
    }
  }

  /* Calculate conditional node probabilities*/
  *p_nAB = (double)(*ncom_nod) / (double)(*nB_nod);
  *p_nBA = (double)(*ncom_nod) / (double)(*nA_nod);

  /* Calculate number of links (that could exist in both networks,
     that is, disregarding those whose endpoint nodes are not in both
     networks) and link overlap */
  pA = netA;
  while ((pA = pA->next) != NULL) {
    lA = pA->neig;
    while ((lA = lA->next) != NULL) {
      /* If both nodes exist in both networks... */
      if (((n1 = GetNodeDict(pA->label, node_dictB)) != NULL) &&
	  ((n2 = GetNodeDict(lA->ref->label, node_dictB)) != NULL)) {
	ntot_lin++;
	(*nA_lin) += 1;
	if (IsThereLink(n1, n2) == 1) {
	  (*ncom_lin) += 1;
	}
	else {
	  fprintf(stderr, ">> %s %s\n", n1->label, n2->label);
	}
      }
    }
  }

  pB = netB;
  while ((pB = pB->next) != NULL) {
    lB = pB->neig;
    while ((lB = lB->next) != NULL) {
      /* If both nodes exist in both networks... */
      if (((n1 = GetNodeDict(pB->label, node_dictA)) != NULL) &&
	  ((n2 = GetNodeDict(lB->ref->label, node_dictA)) != NULL)) {
	(*nB_lin) += 1;
	if (IsThereLink(n1, n2) != 1) {
	  ntot_lin++;
	  fprintf(stderr, "<< %s %s\n", n1->label, n2->label);
	}
      }
    }
  }

  (*nA_lin) /= 2;
  (*nB_lin) /= 2;
  (*ncom_lin) /= 2;
      
  /* Calculate node conditional probabilities */
  *p_lAB = (double)(*ncom_lin) / (double)(*nB_lin);
  *p_lBA = (double)(*ncom_lin) / (double)(*nA_lin);

  /* Free memory and return */
  FreeLabelDict(node_dictA);
  FreeLabelDict(node_dictB);
  return;
}

/*
  -----------------------------------------------------------------------------
  Given two networks, the subroutine returns the "sum" of the
  networks. The sum network contains all nodes and links that appear
  in netA and/or netB.
  -----------------------------------------------------------------------------
*/
struct node_gra *
AddTwoNetworks(struct node_gra *netA, struct node_gra *netB)
{
  struct node_gra *netSum = CopyNetwork(netA);
  void *nodeDictA, *nodeDictSum;
  struct node_gra *p;
  struct node_lis *n;

  /* Add nodes in B missing from A */
  nodeDictA = MakeLabelDict(netA);
  p = netB;
  while((p = p->next) != NULL)
    if (GetNodeDict(p->label, nodeDictA) == NULL)
      CreateNodeGraph(netSum, p->label);
  FreeLabelDict(nodeDictA);

  /* Add all links in B (repetitions are taken care of by
     AddAdjacency) */
  nodeDictSum = MakeLabelDict(netSum);
  p = netB;
  while ((p = p->next) != NULL) {
    n = p->neig;
    while ((n = n->next) != NULL) {
      AddAdjacency(GetNodeDict(p->label, nodeDictSum),
		   GetNodeDict(n->ref->label, nodeDictSum),
		   1, 0, 1, 0);
    }
  }
  FreeLabelDict(nodeDictSum);

  /* Done */
  return netSum;
}
