#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <search.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "partition.h"
#include "io.h"

int
node_compare(const void *node1, const void *node2)
{
    return strcmp(((const struct node *) node1)->label,
        ((const struct node *) node2)->label);
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
@brief Free the memory allocated to binary tree.
*/
void
free_node_tree(void *root)
{
  struct node_t *focal = (struct node_t*) root;
  if (focal != NULL){
	struct node *leaf = (struct node*) focal->key;
   	free_node_tree(focal->left);
   	free_node_tree(focal->right);
	free(leaf);
	free(focal);
  }
}


/**
Read an edge list from a file and convert it to arrays.
@param nodes_in_p,nodes_out_p,weights_p Edge list, nodes id and weight of each edge.
@param labels_p Labels of nodes whose index are in nodes_in/out
@param E_p N_p Number of edges, number of nodes.
@param weighted If True
@param bipartite 0,1,2. If 1 (resp 2), treat the conlumns independently and store only the first (resp 2nd) column labels.
**/
int
EdgeListFileInput(FILE *inFile, int weighted, int bipartite, unsigned int **nodes_in_p,
                  unsigned int **nodes_out_p,double **weights_p,
							    char ***labels_p, unsigned int *E_p, unsigned int *N_p){

  unsigned int *nodes_in=NULL, *nodes_out=NULL;
  double *weights=NULL;
  char **labels=NULL;

  if (bipartite > 2) return 2;

  int i = 0;

  int noReadItems;
  char *line = NULL;
  size_t bufsiz = 0;
  ssize_t nbytes;

  unsigned int Nmax = 10;
  unsigned int Emax = Nmax*Nmax;

  labels = calloc(Nmax+1,sizeof(char*));
  weights = malloc((Emax+1)*sizeof(double));
  nodes_in = malloc((Emax+1)*sizeof(unsigned int));
  nodes_out = malloc((Emax+1)*sizeof(unsigned int));

  unsigned int N = 0;
  unsigned int N_other_side = 0;
  unsigned int E = 0;

  struct node *nodes_root = NULL;
  struct node *nodes_root_bipartite = NULL;

  struct node *nodeptr = NULL;

  unsigned int noReadItemsExpected = weighted ? 3:2;

  // Go through the input file
  while ((nbytes = getline(&line, &bufsiz, inFile)) != -1){
    //printf("N=%d/%d, E=%d/%d\n",N,Nmax,E,Emax);

    struct node *node1 = malloc(sizeof(struct node));
    struct node *node2 = malloc(sizeof(struct node));
    struct node *temp = NULL;


    /* Read the labels (and edge weight, if necessary) */
    if (weighted == 0) {
      noReadItems = sscanf(line, "%s %s", node1->label, node2->label);
      weights[E] = 1;
	  } else {
      noReadItems = sscanf(line,"%s %s %lf", node1->label, node2->label, &weights[E]);
    }

    if (bipartite>1){
      temp = node1;
      node1 = node2;
      node2 = temp;
    }

    /* Check input sanity */
    if(noReadItems != noReadItemsExpected){
      printf ("Failed to read input: not enough fields in line %s (%d!=%d). \n",
              line, noReadItems,noReadItemsExpected);
      return 1;
    }

    // Check if the first label was encountered already.
    nodeptr = (struct node*) tsearch((void*) node1, (void **)&nodes_root, node_compare);
    struct node *answer = NULL;
    answer = *(struct node **)nodeptr;

    // If the node was not present in the tree, save a new label,
    // Otherwise free the read node to avoid memory leak.
    if (answer==node1){
      answer->id = N;
      labels[N] = malloc(sizeof(answer->label));
      strcpy(labels[N],answer->label);
      N++;
    }else{
      free(node1);
    }
    nodes_in[E] = answer->id;

    // Likewise for the second label
    struct node *answer2 = NULL;
    if(!bipartite){
      nodeptr = (struct node*) tsearch((void*) node2, (void **)&nodes_root, node_compare);
      answer2 = *(struct node **)nodeptr;
      if (answer2==node2){
        answer2->id = N;
        labels[N] = malloc(sizeof(answer2->label));
        strcpy(labels[N],answer2->label);
        N++;
      }else{
        free(node2);
      }
    }
    // For bipartite networks,
    // - Each colums is a different namespace.
    // - the node label of the second column are not saved.
    else{
      nodeptr = (struct node*) tsearch((void*) node2, (void **)&nodes_root_bipartite, node_compare);
      answer2 = *(struct node **)nodeptr;
      if (answer2==node2){
        answer2->id = N_other_side;
        N_other_side++;
      } else{
        free(node2);
      }
    }
    nodes_out[E] = answer2->id;
    E++;
    // Increase the allocated memory if needed.
    if(N+2>=Nmax){
      printf("%d > %d, doubling memory for nodes\n",N,Nmax );
      Nmax *= 2;
      labels = realloc(labels, (Nmax+1)*sizeof(char*));
    }
    if(E+1>=Emax){
      printf("%d > %d, doubling memory for edges",E,Emax );
      Emax *= 2;
      nodes_in = realloc(nodes_in, (Emax+1)*sizeof(unsigned int));
      nodes_out = realloc(nodes_out, (Emax+1)*sizeof(unsigned int));
      weights = realloc(weights, (Emax+1)*sizeof(double));
    }
  }

  // Free unused memory
  nodes_in = realloc(nodes_in, (E+1)*sizeof(unsigned int));
  nodes_out = realloc(nodes_out, (E+1)*sizeof(unsigned int));
  weights = realloc(weights, (E+1)*sizeof(double));
  labels = realloc(labels, (N+1)*sizeof(char*));
  free_node_tree(nodes_root);
  free_node_tree(nodes_root_bipartite);

  // Assign
  *nodes_in_p = nodes_in;
  *nodes_out_p = nodes_out;
  *labels_p = labels;
  *weights_p = weights;
  *E_p = E;
  *N_p = N;
  free(line);

  if (!bipartite)
    printf("Read %d nodes and %d edges \n ", *N_p,*E_p);
  else
    printf("Read %d/%d nodes and %d edges \n ", *N_p,N_other_side,*E_p);
  //
  // for (i=0;i<*E_p;i++){
  //   printf("%d\t [%d] %s \t [%d] %s \t %f\n",i, nodes_in[i],labels[nodes_in[i]],nodes_out[i],labels[nodes_out[i]],weights[i]);
  // }
  // for (i=0;i<*N_p+1;i++){
  //   printf("%d %s\n",i, labels[i]);
  // }
  return 0;
}



void
TabularOutput(FILE *outf,
			  char **labels,
			  Partition *part,
			  double *connectivity,
			  double *participation){
  unsigned int i=0;
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
  unsigned int mod;
  Node *node = NULL;
  for(mod=0; mod<part->M; mod++){
	for(node = part->modules[mod]->first; node!=NULL; node=node->next){
	  fprintf(outf,"%s\t",labels[node->id]);
	}
	fprintf(outf,"\n");
  }
}

int
AssignNodesToModulesFromFile(FILE *inF,
							 Partition *part,
							 char **labels){
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
	  i = (int) ep->data;
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
