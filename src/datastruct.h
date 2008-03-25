/*
  datastruct.h
  $LastChangedDate$
  $Revision$
*/

#include <stdbool.h>

#ifndef RGRAPH_DATASTRUCT_H
#define RGRAPH_DATASTRUCT_H 1

/*
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
  STACK
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
*/
/*
  -----------------------------------------------------------------------------
  Definition of the stack structure
  -----------------------------------------------------------------------------
*/
struct stack{
  struct stack_node *top;   // first element in the stack
  int length;               // number of elements in the stack
};

/*
  -----------------------------------------------------------------------------
  Definition of the node_stack structure
  -----------------------------------------------------------------------------
*/
struct stack_node{
  struct stack_node *next;   // next element in the stack
  void *data;                // the data being stored in the node
};

/*
  -----------------------------------------------------------------------------
  Stack functions
  -----------------------------------------------------------------------------
*/
struct stack *stack_create();
bool stack_empty(struct stack *theStack);
void stack_push(struct stack *theStack, void *theElement);
void *stack_pop(struct stack *theStack);
void *stack_top(struct stack *theStack);
void stack_clear(struct stack *theStack);
void stack_free(struct stack *theStack);

/*
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
  QUEUE
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
*/
/*
  -----------------------------------------------------------------------------
  Definition of the queue structure
  -----------------------------------------------------------------------------
*/
struct queue{
  struct queue_node *first;   // first element in the queue
  struct queue_node *last;    // last element in the queue
  int length;                 // number of elements in the queue
};

/*
  -----------------------------------------------------------------------------
  Definition of the node_queue structure
  -----------------------------------------------------------------------------
*/
struct queue_node{
  struct queue_node *next;   // next element in the queue
  void *data;                // the data being stored in the node
};

/*
  -----------------------------------------------------------------------------
  Queue functions
  -----------------------------------------------------------------------------
*/
struct queue *queue_create();
bool queue_empty(struct queue *theQueue);
void queue_enqueue(struct queue *theQueue, void *theElement);
void *queue_dequeue(struct queue *theQueue);
void *queue_first(struct queue *theQueue);
void queue_clear(struct queue *theQueue);
void queue_free(struct queue *theQueue);


#endif /* !RGRAPH_DATASTRUCT_H */
