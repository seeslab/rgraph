/*
  datastruct.c
  $LastChangedDate$
  $Revision$
*/

#include <stdlib.h>

#include "datastruct.h"

/*
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
  Stack functions
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
*/
/*
  -----------------------------------------------------------------------------
  Create a stack
  -----------------------------------------------------------------------------
*/
struct stack *
stack_create()
{
  struct stack *theStack;
  theStack = (struct stack *)calloc(1, sizeof(struct stack));
  theStack->top = NULL;
  theStack->length = 0;

  return theStack;
}

/*
  -----------------------------------------------------------------------------
  Return true if stack is empty
  -----------------------------------------------------------------------------
*/
bool
stack_empty(struct stack *theStack)
{
  if (theStack->top == NULL)
    return true;
  else
    return false;
}

/*
  -----------------------------------------------------------------------------
  Push element onto stack
  -----------------------------------------------------------------------------
*/
void
stack_push(struct stack *theStack, void *theElement)
{
  struct stack_node *theStackNode;

  theStackNode = (struct stack_node *)calloc(1, sizeof(struct stack_node));
  theStackNode->data = theElement;
  theStackNode->next = theStack->top;
  theStack->top = theStackNode;
  theStack->length += 1;

  return;
}

/*
  -----------------------------------------------------------------------------
  Return top node and remove it from the stack
  -----------------------------------------------------------------------------
*/
void *
stack_pop(struct stack *theStack)
{
  void *theData;
  struct stack_node *temp;

  if (stack_empty(theStack)) {
    return NULL;
  }
  else {
    temp = theStack->top;
    theData = temp->data;
    theStack->top = temp->next;
    free(temp);
    theStack->length -= 1;
    return theData;
  }
}

/*
  -----------------------------------------------------------------------------
  Return top node
  -----------------------------------------------------------------------------
*/
void *
stack_top(struct stack *theStack)
{
  if (stack_empty(theStack))
    return NULL;
  else
    return theStack->top->data;
}

/*
  -----------------------------------------------------------------------------
  Clear stack
  -----------------------------------------------------------------------------
*/
void
stack_clear(struct stack *theStack)
{
  while (!stack_empty(theStack))
    stack_pop(theStack);
  return;
}

/*
  -----------------------------------------------------------------------------
  Free memory allocated for a stack
  -----------------------------------------------------------------------------
*/
void
stack_free(struct stack *theStack)
{
  stack_clear(theStack);
  free(theStack);
  return;
}

/*
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
  Queue functions
  -----------------------------------------------------------------------------
  -----------------------------------------------------------------------------
*/
/*
  -----------------------------------------------------------------------------
  Create a queue
  -----------------------------------------------------------------------------
*/
struct queue *
queue_create()
{
  struct queue *theQueue;
  theQueue = (struct queue *)calloc(1, sizeof(struct queue));
  theQueue->first = NULL;
  theQueue->last = NULL;
  theQueue->length = 0;

  return theQueue;
}

/*
  -----------------------------------------------------------------------------
  Return true if queue is empty
  -----------------------------------------------------------------------------
*/
bool
queue_empty(struct queue *theQueue)
{
  if (theQueue->first == NULL)
    return true;
  else
    return false;
}

/*
  -----------------------------------------------------------------------------
  Push element onto queue
  -----------------------------------------------------------------------------
*/
void
queue_enqueue(struct queue *theQueue, void *theElement)
{
  struct queue_node *theQueueNode;

  theQueueNode = (struct queue_node *)calloc(1, sizeof(struct queue_node));
  theQueueNode->data = theElement;
  theQueueNode->next = NULL;
  
  if (queue_empty(theQueue)) {
    theQueue->first = theQueue->last = theQueueNode;
  }
  else {
    theQueue->last->next = theQueueNode;
    theQueue->last = theQueueNode;
  }

  theQueue->length += 1;

  return;
}

/*
  -----------------------------------------------------------------------------
  Return first node and remove it from the queue
  -----------------------------------------------------------------------------
*/
void *
queue_dequeue(struct queue *theQueue)
{
  void *theData;
  struct queue_node *temp;

  if (queue_empty(theQueue)) {
    return NULL;
  }
  else {
    temp = theQueue->first;
    theData = temp->data;
    theQueue->first = temp->next;
    free(temp);
    theQueue->length -= 1;
    return theData;
  }
}

/*
  -----------------------------------------------------------------------------
  Return first node
  -----------------------------------------------------------------------------
*/
void *
queue_first(struct queue *theQueue)
{
  if (queue_empty(theQueue))
    return NULL;
  else
    return theQueue->first->data;
}

/*
  -----------------------------------------------------------------------------
  Clear queue
  -----------------------------------------------------------------------------
*/
void
queue_clear(struct queue *theQueue)
{
  while (!queue_empty(theQueue))
    queue_dequeue(theQueue);
  return;
}

/*
  -----------------------------------------------------------------------------
  Free memory allocated for a queue
  -----------------------------------------------------------------------------
*/
void
queue_free(struct queue *theQueue)
{
  queue_clear(theQueue);
  free(theQueue);
  return;
}
