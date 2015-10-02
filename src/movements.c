#include "partition.h"
#include "movements.h"


/**
The variation in modularity induced by moving a node i from a group G
to a group H is given by:

2 * (  k_{i,H} - k_{i,G} + k_i (K_G - K_H) )

with
- k_i the strength of i,
- k_(i,H) the sum of the strength from i to H,
- k_H the sum of the strength of nodes of H.

@param nodeid the node to move.
@param newModuleid the target module.
**/
double dEChangeModule(unsigned int nodeid,
					  unsigned int newModuleid,
					  Partition *part, AdjaArray *adj){
  double dE = 0.0;
  unsigned int i = 0;
  unsigned int old = part->nodes[nodeid]->module;
  Module * OldModule = part->modules[old];
  Module * NewModule = part->modules[newModuleid];
  Node * Node = part->nodes[nodeid];

  // Loop through the neighbors and compute the strength to
  // old and new group.
  for (i=adj->idx[nodeid]; i<adj->idx[nodeid+1]; i++){
	int j = adj->neighbors[i]; // The index of the neighbor.
	if (part->nodes[j]->module == old)
	  dE -= adj->strength[i]; // - k_(i,G)
	else if (part->nodes[j]->module == newModuleid)
	  dE += adj->strength[i]; // + k_(i,H)
  }

  // Group properties: k_i * ( K_H-k_i - K_G)
  dE += Node->strength * (OldModule->strength - Node->strength  - NewModule->strength);
  return(2*dE);
}

/**
The variation in modularity induced by merging modules G and H is
given by:

2 * (  SUM_FOR_ALL_i_IN_G[k_{i,H}] - (K_G * K_H) )

with
- k_i the strength of i,
- k_(i,H) the sum of the strength from i to H,
- k_H the sum of the strength of nodes of H.
**/
double
dEMergeModules(unsigned int moduleId1,
			   unsigned int moduleId2,
			   Partition *part, AdjaArray *adj){

  Module *small;
  unsigned int idlarge;
  Node *node, **nodes;
  int i;
  double dE=0;
  nodes = part->nodes;

  // Get the smaller and larger module.
  if (part->modules[moduleId1]->size
	  > part->modules[moduleId2]->size){
	small = part->modules[moduleId2];
	idlarge = moduleId1;
  }else{
	small = part->modules[moduleId1];
	idlarge = moduleId2;
  }

  // Compute the strength to the large module of all nodes in the
  // small module... (SUM_FOR_ALL_i_IN_G[k_{i,H}])
  for(node=small->first; node!=NULL; node = node->next){
	// Loop through neigbors
	for (i=adj->idx[node->id]; i<adj->idx[(node->id)+1]; i++){
	  int j = adj->neighbors[i]; // The index of the neighbor.
	  if (nodes[j]->module == idlarge )
		dE += adj->strength[i];
	}
  }

  // Product of the module strength. (-K_G*K_H)
  dE -= (part->modules[moduleId1]->strength
		 * part->modules[moduleId2]->strength);
  return(2*dE);
}


/**
 Change the module of a node within a partition.

This involve, updating its properties, bookkeeping the module and
partition properties and updating pointers in the doubly linked-list.
**/
void
ChangeModule(unsigned int nodeId,
			 unsigned int newModuleId,
			 Partition *part){
  Node *node = part->nodes[nodeId];
  unsigned int oldModuleId = part->nodes[nodeId]->module;
  Module *OldModule = part->modules[oldModuleId];
  Module *NewModule = part->modules[newModuleId];

  // Changing the node property.
  node->module = NewModule->id;

  // Bookkeeping modules and partition properties.
  OldModule->strength -= node->strength;
  OldModule->size--;
  if (!OldModule->size){
	part->nempty++;
	//OldModule->strength=0;
  }
  if (!NewModule->size){
	part->nempty--;
  }
  NewModule->strength += node->strength;
  NewModule->size++;

  // Moving the node in the double linked list.

  // Pr <-> N <-> Nx => Pr <-> Nx.
  // Extract the node from its old doubly linked list.
  if (node->prev != NULL)
	node->prev->next = node->next;
  else
	OldModule->first = node->next;
  if (node->next != NULL)
	node->next->prev = node->prev;
  else
	OldModule->last = node->prev;

  // G <-> N' => G <-> N <-> N'
  node->next = NewModule->first;
  node->prev = NULL;
  if (NewModule->first != NULL)
	NewModule->first->prev = node; //Used to be NULL;
  else
	NewModule->last = node;
  NewModule->first = node;
}


/**
 Merge two modules of a partition.


- Whithout effect if one of the modules is empty. Otherwise, append
  the smaller to the begining of the bigger.

@param id1,id2 Modules to merge.
@param part Partition containing those modules.
**/
void
MergeModules(unsigned int id1, unsigned int id2,
			 Partition *part){

  Module *small, *large;
  Node *node;

  // If one of the modules is empty, stop here.
  if (!part->modules[id1]->size || !part->modules[id2]->size)
	return;

  // Get the smaller and larger module. Because we want to append the
  // smaller to the begining of the longer linked list.
  if (part->modules[id1]->size > part->modules[id2]->size){
	small = part->modules[id2];
	large = part->modules[id1];
  }else{
	small = part->modules[id1];
	large = part->modules[id2];
  }
  // Updating nodes properties...
  for(node=small->first; node!=NULL; node = node->next)
	node->module = large->id;

  // Bookkeeping partition properties.
  part->nempty++;

  // Bookkeeping modules properties.
  large->size += small->size;
  small->size = 0;
  large->strength += small->strength;
  small->strength = 0;

  // Moving the linked lists.
  large->first->prev = small->last;
  small->last->next = large->first;
  large->first = small->first;
  small->last = NULL;
  small->first = NULL;
}

/**
Try to split a target module in two along the connected component in
the induced subgraph. This function return the number of connected
components.

If the subgraph spawning from targetModuleId has at least two
connected components, this function guarantee that at least one of
them will be moved into emptyModuleId. Additional connected components
are randomly assigned to one of the two modules. If there is only one
connected component, the split is trivial.

@param targetModuleId The module to split.
@param emptyModuleId The module to split into.
@param part The partition.
@param adj The adjacency array to traverse to look for connected components.
@param gen Random number generator.
@return The number of connected components. If it is 1, all nodes are
		still in the targetModule, and this function was without effect.
**/
unsigned int
SplitModuleByComponent(unsigned int targetModuleId,
					   unsigned int emptyModuleId,
					   Partition *part,
					   AdjaArray *adj,
					   gsl_rng *gen){
  Module *target = part->modules[targetModuleId];
  Stack *to_visit, *to_move;
  Node *node;
  unsigned int to_find = target->size;
  unsigned int components = 0;
  unsigned int add_to_empty, j, i, k;
  unsigned int *visited = calloc(part->N,sizeof(unsigned int));
  if (visited==NULL)
	  perror("Error while splitting module by connected component");

  to_visit = CreateStack(part->N);
  to_move = CreateStack(target->size);
  to_find = target->size;

  // Loop trough the nodes in the module.
  for(node=target->first; node!=NULL && to_find; node = node->next){
	if(!visited[node->id]){

	  visited[node->id]=1; //Now the node is visited

	  // If it was an unvisited node, it means we are starting a new
	  // connected component.
	  components++;
	  if(components==1) //If it is the first con. compo., we keep it in its module.
		add_to_empty = 0;
	  else if(components==2) //If it is the second we switch to the other module.
		add_to_empty = 1;
	  else //Otherwise, we pick a module at random.
		add_to_empty = gsl_rng_uniform(gen)>.5? 1 : 0;

	  // ADD the node to the queue of nodes to visit.
	  AddToStack(node->id,to_visit);

	  // Traverse the graph.
	  while((j=PopFromStack(to_visit))!=-1 && to_find){
		// Is it one of our target nodes ?
		if (part->nodes[j]->module == targetModuleId){
		  to_find--;
		  //If we are supposed to move it, add it to the list of nodes
		  //to move.
		  if(add_to_empty)
			AddToStack(j,to_move);
		}

		//loop through neighbors and add them to the queue.
		for (i=adj->idx[j];i<adj->idx[j+1];i++){
		  k = adj->neighbors[i]; // The index of the neighbor.

		  if(!visited[k]){
			visited[k] = 1;
			// Add the node to the queue of nodes to visit.
			AddToStack(k,to_visit);
		  }
		} // End of loop through neighbors.

	  } // END of loop through queue.
	} // End of if node in module unvisited.
  } // End of loop through module.

  // Move the nodes we marked for moving.
  while((i=PopFromStack(to_move))!=-1)
	ChangeModule(i,emptyModuleId,part);

  FreeStack(to_move);
  FreeStack(to_visit);
  free(visited);
  return(components);
}

/**
Allocate the memory for a stack of maximum size N.
 **/
Stack *
CreateStack(unsigned int N){
  Stack * st = (Stack *) malloc(sizeof(Stack));
  st->maxsize = N;
  st->top = -1;
  st->items = (int*) calloc(N,sizeof(int));
  if (st==NULL || st->items==NULL)
	  perror("Error while allocating stack");
  return(st);

}

/**
Free the memory allocated by a stack
**/
void
FreeStack(Stack *st){
  free(st->items);
  free(st);
  st = NULL;
}

/**
Add a value to the the top of stack.
 **/
void
AddToStack(unsigned int value, Stack * st){
  st->top++;
  st->items[st->top] = value;
}

/**
Get the value at the top of the stack. Returns -1 if the stack is
empty.
 **/
int
PopFromStack(Stack *st){
  unsigned int value;
  if(st->top>=0){
	value = st->items[st->top];
	st->top--;
	return(value);
  }
  else
	return(-1);
}
