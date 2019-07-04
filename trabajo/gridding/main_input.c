#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include "hazmath.h"
/**********************************************************************/
void input_grid_free(input_grid *g);
void input_grid_print(input_grid *g);
void coo2csr(INT nrow,INT ncol,INT nnz,					\
	     INT *row_idx,INT *col_idx, void *aval,			\
	     INT *ia,INT *ja, void *bval,				\
	     size_t elsize);
char **splits(char *s, const char *d, INT *num);
void read_data(char *data_coordsystems,		\
	       char *data_vertices,		\
	       char *data_edges,		\
	       input_grid *g);
void get_out(char *pattern, size_t le);
char *make_string_from_file(FILE *the_file, size_t *length_string);
char *get_substring(char *pattern,		\
		    size_t *length_substring,	\
		    char *the_string);
input_grid *parse_input_grid(const char *input_file_grid);
/***************************************************************/
void set_input_grid(input_grid *g)
{
  INT i,j,iri,ici;
  INT *p=calloc(2*g->nv,sizeof(INT));// permutation and inverse permutation;
  dlexsort(g->nv, g->dim,g->x,p);
  //  for (i=0;i<g->nv;i++)fprintf(stdout,"\n%d-->%d",p[i],i);  
  //  fprintf(stdout,"\n"); 
  // use p as a working array to store labels;
  INT *invp=p+g->nv; // inverse permutation;
  for (i=0;i<g->nv;i++){
    invp[p[i]]=i;
    p[i]=g->labels[i];
  }
  /* permute labels (coordinate systems for vertices */
  for (i=0;i<g->nv;i++)
    g->labels[i]=p[invp[i]]; // fix coord systems;
  for (i=0;i<g->ne;i++){
    iri=invp[g->seg[3*i]];
    ici=invp[g->seg[3*i+1]];
    //    fprintf(stdout,"\n(%d,%d)-->[%d,%d]: div=%d",g->seg[3*i],g->seg[3*i+1],iri,ici,g->seg[3*i+2]);
    if(iri<ici){
      g->seg[3*i]=iri;
      g->seg[3*i+1]=ici;
    } else {
      g->seg[3*i]=ici;
      g->seg[3*i+1]=iri;
    }
    /* set up divisions */
    j=g->seg[3*i+1]-g->seg[3*i]; // should be always positive;
    if(g->seg[3*i+2]>p[j])
      p[j]=g->seg[3*i+2];
  }  
  for (i=0;i<g->ne;i++){
    j=g->seg[3*i+1]-g->seg[3*i];
    g->seg[3*i+2]=p[j]; 
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  if(p) free(p); // this also frees invp;
  return;
}
/*********************************************************************/
int main(int argc, char **argv){
  char input_grid_file[256]={"grid.input"};
  input_grid *g=parse_input_grid(input_grid_file);
  set_input_grid(g);
  //
  input_grid_print(g);
  input_grid_free(g);
  return 0;
}
