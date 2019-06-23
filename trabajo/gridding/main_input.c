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
/* void input_grid_free(input_grid *g); */
/* void input_grid_print(input_grid *g); */
/* void coo2csr(INT nrow,INT ncol,INT nnz,					\ */
/* 	     INT *row_idx,INT *col_idx, void *aval,			\ */
/* 	     INT *ia,INT *ja, void *bval,				\ */
/* 	     size_t elsize); */
/* char **splits(char *s, const char *d, INT *num); */
/* void read_data(char *data_coordsystems,		\ */
/* 	       char *data_vertices,		\ */
/* 	       char *data_edges,		\ */
/* 	       input_grid *g); */
/* void get_out(char *pattern, size_t le); */
/* char *make_string_from_file(FILE *the_file, size_t *length_string); */
/* char *get_substring(char *pattern,		\ */
/* 		    size_t *length_substring,	\ */
/* 		    char *the_string); */
/* input_grid *parse_input_grid(const char *input_file_grid); */
/********************************************************************/
void lexsort(const INT nr, const INT nc,REAL *a,INT *p);
/***************************************************************/
void set_input_grid(input_grid *g)
{
  INT i,j,k,iri,ici;
  INT *p=calloc(2*g->nv,sizeof(INT));// permutation and inverse permutation;
  lexsort(g->nv, g->dim,g->x,p);
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
  for (i=0;i<g->seg->nnz;i++){
    iri=invp[g->seg->rowind[i]];
    ici=invp[g->seg->colind[i]];
    if(iri<ici){
      g->seg->rowind[i]=iri;
      g->seg->colind[i]=ici;
    } else {
      g->seg->rowind[i]=ici;
      g->seg->colind[i]=iri;
    }
    //    fprintf(stdout,"\n(%d,%d)-->[%d,%d]: div=%d",g->seg->rowind[i],g->seg->colind[i],iri,ici,g->seg->val[i]);
    /* set up divisions */
    j=g->seg->colind[i]-g->seg->rowind[i]; // should be always positive;
    if(g->seg->val[i]>p[j])
      p[j]=g->seg->val[i];
  }  
  for (i=0;i<g->seg->nnz;i++){
    j=g->seg->colind[i]-g->seg->rowind[i];
    g->seg->val[i]=p[j]; 
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg->rowind[i],g->seg->colind[i],g->seg->val[i]);
  }
  if(p) free(p); // this also frees invp;
  return;
}
/*********************************************************************/
int main(int argc, char **argv){
  char input_grid_file[256]={"grid.input"};
  input_grid *g=parse_input_grid(input_grid_file);
  input_grid_print(g);
  set_input_grid(g);
  input_grid_print(g);
  input_grid_free(g);
  return 0;
}
