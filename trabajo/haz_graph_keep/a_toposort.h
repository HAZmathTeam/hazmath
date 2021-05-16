//////////////////////////////////////////////////////////////
// A structure to represent a stack
// adjacency graph
typedef struct graphij{    
  int nv;
  int ne;    
  int *i;
  int *j;
  double *weights;    // we should make this a matrix as well.
  double *bweights; 
} graphij; 

typedef struct grapha{    
  int nv;
  int ne;
  int *ia;
  int *ja;
  double *weights;    // we should make this a matrix
  double *bweights;    
} grapha;

typedef struct stack {
  int top;
  int num_items;
  int null_item;
  int *items;
} stack;
//////////////////////////////////////////////////////////////
static void print_graphij(graphij *g){
  int i;
  if(g->nv<=0 || g->ne<=0) return;
  fprintf(stdout,"\n===============\n");
  fprintf(stdout,"Printing coo graph: %s\n",__FUNCTION__);
  fprintf(stdout,"\nNV=%d, NE=%d",g->nv,g->ne);
  for (i=0;i<g->ne;i++){
    fprintf(stdout,"\n[%d-->%d]: weight=%5.2f; bweight=%5.2f",	\
    	    g->i[i],g->j[i],g->weights[i],g->bweights[i]);
  }
  fprintf(stdout,"\n");
  fprintf(stdout,"\n===============\n");
  fflush(stdout);
  return;
}
static void print_grapha(grapha *g){
  int i,ij,iaa,iab;
  if(g->nv<=0 || g->ne<=0) return;
  fprintf(stdout,"\n===============\n");
  fprintf(stdout,"Printing adj graph: %s\n",__FUNCTION__);
  fprintf(stdout,"\nNV=%d, NE=%d",g->nv,g->ne);
  for (i=0;i<g->nv;i++){
    iaa=g->ia[i];iab=g->ia[i+1];
    fprintf(stdout,"\nrow    %5d: ",i);
    for (ij=iaa;ij<iab;ij++){
      fprintf(stdout," %5d",g->ja[ij]);
    }
    fprintf(stdout,"\nvalues      : ");
    for (ij=iaa;ij<iab;ij++){
      fprintf(stdout," %5.2f(%5.2f)",g->weights[ij],g->bweights[i]);
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n===============\n");
  fflush(stdout);
  return;
}
/******************************************************************************************/
graphij *graphij_init(const int nv)
{
  graphij *g=malloc(sizeof(graphij));
  g->nv=nv;
  g->ne=0;    
  g->i=NULL;
  g->j=NULL;
  g->weights=NULL;
  g->bweights=NULL;
  return g;
}
/******************************************************************************************/
void edge_add(graphij *g, const int a, const int b, const double w, const double bw)
{  
  int ne=g->ne;
  double tmp;
  g->ne++;    
  g->i=realloc(g->i,g->ne*sizeof(int));
  g->j=realloc(g->j,g->ne*sizeof(int));
  g->weights=realloc(g->weights,g->ne*sizeof(double));
  g->bweights=realloc(g->bweights,g->ne*sizeof(double));
  g->i[ne]=a;
  g->j[ne]=b;
  g->weights[ne]=w;
  g->bweights[ne]=bw;
}
/******************************************************************************************/
void coo2adj ( graphij *g, grapha *adj)
// from hazmath: converts (i,j,value) to adjacency (CSR) format.
{   
    // get size
    const int nv=g->nv, ne=g->ne;
    adj->ia=(int *)calloc(nv+1,sizeof(int));
    adj->ja=(int *)calloc(ne,sizeof(int));   
    adj->weights=(double *)calloc(ne,sizeof(double));
    adj->bweights=(double *)calloc(ne,sizeof(double));
    adj->nv=nv;adj->ne=ne;
    // local variables
    int *ia = adj->ia;
    int *ja = adj->ja;
    double *adjval = adj->weights;
    double *adjbval = adj->bweights;
    int *row_idx = g->i;
    int *col_idx = g->j;
    double *gval = g->weights;
    double *gbval = g->bweights;
    int i, iind, jind;    
    int *ind = (int *) calloc(nv+1,sizeof(int));    
    // initialize
    for(i=0;i<=nv;++i) ind[i]=0;
    //    memset(ind, 0, sizeof(int)*(nv+1));   
    // count number of nonzeros in each row
    for (i=0; i<g->ne; ++i) ind[row_idx[i]+1]++;    
    // set row pointer
    ia[0] = 0;
    for (i=1; i<=nv; ++i) {
        ia[i] = ia[i-1]+ind[i];
        ind[i] = ia[i];
    }    
    // set column index and values
    for (i=0; i<g->ne; ++i) {
        iind = row_idx[i];
        jind = ind[iind];
        ja[jind] = col_idx[i];
        adjval[jind] = gval[i];
        adjbval[jind] = gbval[i];
        ind[iind] = ++jind;
    }    
    if (ind) free(ind);
    return;
}
/*****************************************************************************************/
// function to create a stack of given num_items. It initializes size of
// stack as 0
stack* stack_init(int num_items)
{
  stack *s = (stack *)malloc(sizeof(stack));
  s->num_items = num_items;
  s->top = -1;
  s->null_item=INT_MIN;
  s->items = (int*)malloc(s->num_items * sizeof(int));
  return s;
}
// Function to push an item to stack. It increases top by 1
void push_in(stack* s, int item)
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
int pop_out(stack* s)
{
  if (s->top==-1)
    return s->null_item;
  int top_el=s->items[s->top];
  s->top--;
  return top_el;
}
// Function to return the top from stack without removing it
int get_top(stack *s)
{
  if (s->top==-1) return s->null_item;
  return s->items[s->top];
}
// A recursive function used by topo_sort
void topo_sort_one(int v, int *ia, int *ja, short *mask,stack *s)
{
  int i,vrow;
  // Mark the current node as visited.
  mask[v] = (short )1;
  for (vrow=ia[v];vrow<ia[v+1]; ++vrow){
    i=ja[vrow];
    if (!mask[i])    
      topo_sort_one(i,ia,ja,mask, s);
  }
    push_in(s,v);
}
// The function to do Topological Sort.
// It uses recursive topo_sort_one()
void topo_sort(grapha *g)
{
  int i,item;  
  stack *s = stack_init(g->nv);
  // Mark all the vertices as not visited
  short *mask = (short *)calloc(g->nv,sizeof(short));
  for (i=0;i<g->nv;++i) mask[i]=(short )0;
  //
  for (i=0;i<g->nv;++i)
    if (!mask[i])
      topo_sort_one(i, g->ia,g->ja, mask, s);
  // Print contents of stack
  while (s->top != -1) {
    item=pop_out(s);
    fprintf(stdout,"%d ",item);
  }
}

