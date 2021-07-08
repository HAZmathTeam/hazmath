#include "hazmath.h"
/*****************************************************************************************/
// function to create a stack of given num_items. It initializes size of
// stack as 0
stack* stack_init(INT num_items)
{
  stack *s = (stack *)malloc(sizeof(stack));
  s->num_items = num_items;
  s->top = -1;
  s->null_item=INT_MIN;
  s->items = (INT*)malloc(s->num_items * sizeof(INT));
  return s;
}
// Function to push an item to stack. It increases top by 1
void push_in(stack* s, INT item)
{
  // push ellement in the stack
  if (s->top >= (s->num_items - 1)) {
    fprintf(stderr,"cannot push furter item in stack: max items reached");
    exit(5);
  }
  s->top++;
  s->items[s->top] = item;
  //  printf("%d pushed to stack\n", item);
  return;
}

// Function to remove an item from stack. It decreases top by 1
INT pop_out(stack* s)
{
  if (s->top==-1)
    return s->null_item;
  INT top_el=s->items[s->top];
  s->top--;
  return top_el;
}
// Function to return the top from stack without removing it
INT get_top(stack *s)
{
  if (s->top==-1) return s->null_item;
  return s->items[s->top];
}
// A recursive function used by topo_sort
static void topo_recurrsive(INT v, INT *ia, INT *ja, SHORT *mask,stack *s)
{
  INT i,vrow;
  // Mark the current node as visited.
  mask[v] = (SHORT )1;
  for (vrow=ia[v];vrow<ia[v+1]; ++vrow){
    i=ja[vrow];
    if (!mask[i])    
      topo_recurrsive(i,ia,ja,mask, s);
  }
    push_in(s,v);
}
// The function to do Topological Sort. It uses recursive topo_recurrsive()
void topo_sort(dCSRmat *a)
{ 
  INT i;  
  stack *s = stack_init(a->row);
  // Mark all the vertices as not visited
  SHORT *mask = (SHORT *)calloc(a->row,sizeof(SHORT));
  for (i=0;i<a->row;++i) mask[i]=(SHORT )0;
  for (i=0;i<a->row;++i)
    if (!mask[i])
      topo_recurrsive(i, a->IA,a->JA,mask, s);
  // trace the contents of stack
  /* while (s->top != -1) { */
  /*   item=pop_out(s); */
  /*   //we can print it if we want    fprintf(stdout,"%d ",item); */
  /* } */
  return;
}

