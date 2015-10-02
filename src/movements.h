#ifndef MOVEMENTS_H__
#define MOVEMENTS_H__
#include "partition.h"

// Node movements
void ChangeModule(unsigned int nodeId,
				  unsigned int newModuleId,
				  Partition *part);
void MergeModules(unsigned int id1, unsigned int id2,
				  Partition *part);
unsigned int SplitModuleByComponent(unsigned int targetModuleId,
									unsigned int emptyModuleId,
									Partition *part, AdjaArray *adj,
									gsl_rng *gen);
// Modularity optimisation
double dEChangeModule(unsigned int nodeid,
					  unsigned int newModuleid,
					  Partition *part, AdjaArray *adj);
double dEMergeModules(unsigned int moduleId1,
					  unsigned int moduleId2,
					  Partition *part, AdjaArray *adj);


/**
A stack of integers, used to keep tracks of nodes to visit
when traversing th graph.

- The function AddToStack allows to add an item to the top of the
  stack.
- The function PopFromStack will return the item at the top of the
  stack and remove it, or -1 if the stack is empty.
**/
typedef struct Stack {
  unsigned int maxsize;
  int top;
  int *items;
} Stack;

Stack * CreateStack(unsigned int N);
void FreeStack(Stack *st);
void AddToStack(unsigned int value, Stack * st);
int PopFromStack(Stack *st);

#endif
