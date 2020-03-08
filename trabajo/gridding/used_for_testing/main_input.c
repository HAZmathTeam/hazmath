#include "header_input.h"
/***********************************************************************/
void input_grid_free(input_grid *g);
void input_grid_print(input_grid *g);
input_grid *parse_input_grid(const char *input_file_grid);
const char **input_strings();
/***********************************************************************/
INT *set_input_grid(input_grid *g)
{
  /* 
     Every edge is put into a subset, i.e. two edges (i1,i2) and (j1,j2)
     are considered equivalent iff (i2-i1)=(j2-j1).  The number of
     divisions in an equivalent set of edges is taken to be the
     largest from the equivalence class.  OUTPUT array is a "dim"
     array and for each direction gives the number of partitions.
  */
  INT i,j,k,iri,ici;
  INT *p=calloc(2*g->nv,sizeof(INT));// permutation and inverse permutation;
  //  dlexsort(g->nv, g->dim,g->xv,p);
  //  for (i=0;i<g->nv;i++)fprintf(stdout,"\n%d-->%d",p[i],i);  
  //  fprintf(stdout,"\n"); 
  // use p as a working array to store labels;
  //XXXXXXXXX  INT *invp=p+g->nv; // inverse permutation;
  //  for (i=0;i<g->nv;i++){
  //    invp[p[i]]=i;
  //    p[i]=g->csysv[i];
  //  }
  /* permute labels (coordinate systems for vertices */
  //no  for (i=0;i<g->nv;i++)
  //no    g->csysv[i]=p[invp[i]]; // fix coord systems;
  for (i=0;i<g->ne;i++){
    iri=g->seg[3*i];
    ici=g->seg[3*i+1];
    /* iri=invp[g->seg[3*i]]; */
    /* ici=invp[g->seg[3*i+1]]; */
    /* fprintf(stdout,"\n(%d,%d)-->[%d,%d]: div=%d",g->seg[3*i],g->seg[3*i+1],iri,ici,g->seg[3*i+2]); */
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
  for (i=0;i<g->ne;i++){
    j=g->seg[3*i+1]-g->seg[3*i];
    g->seg[3*i+2]=p[j]; 
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  ilexsort(g->ne, 3,g->seg,p);
  k=0;
  for (i=0;i<g->ne;i++){
    if(g->seg[3*i]) continue;
    j=g->seg[3*i+1]-g->seg[3*i]-1;
    p[k]=g->seg[3*i+2];
    k++;
    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  p=realloc(p,g->dim*sizeof(INT)); // realloc to dimension g->dim
  for (i=0;i<g->dim;i++){
    fprintf(stdout,"\ndirection:%d; div=%d",i,p[i]);
  }
  return p;
}
/*********************************************************************/
int main(int argc, char **argv){
  char input_grid_file[256]={"grid.input"};
  input_grid *g=parse_input_grid(input_grid_file);
  //YES  INT *nd=set_input_grid(g);
  //
  //YES  input_grid_print(g);
  //YES  input_grid_free(g);
  //YES  free(nd);
  return 0;
}
