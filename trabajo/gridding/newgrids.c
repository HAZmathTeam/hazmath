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
#include "grid_defs.h"
#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif
/********************************************************************/
typedef struct /* general cell complex (like hexagonal, etc) */
{
  INT n; /* the dimension of GC */
  INT nv; /* number of 0-dimensional cells (vertices) */
  INT ne; /* number of 1-dimensional cells (edges) */
  INT ncells; /* number of n-dimensional cells or elements */
  INT *nodes; /* array (depending on faces) to hold the neighbors */
  INT *nbr; /* array (depending on faces) to hold the neighbors */
  INT *bndry; /* nv boundary codes for vertices */
  INT *csys; /* nv coord system for a point */
  REAL *x; /*(nv x n) array to hold the coordinates of vertices all in
	     cartesian coordinates*/
  REAL *xe; /*(ne x n) array to hold the coordinates of midpoints of the edges
	     in cartesian coordinates*/
} gcomplex;
/********************************FINCTIONS:*********************/
void cube2simp_free(cube2simp *c2s);
INT reverse(void *arr,INT length, size_t elsize);
cube2simp *cube2simplex(INT dim);
scomplex *umesh(const INT dim, INT *nd, cube2simp *c2s, const INT intype);
void polar2cart(INT dim, REAL *px, REAL *cx);
REAL interp8(cube2simp *c2s, REAL *u, REAL *ue, REAL *xhat);
REAL interp4(cube2simp *c2s, REAL *u, REAL *xhat);
void unirefine(INT *nd,scomplex *sc);
/***************************************************************/
/********************************************************************/
void lexsort(const INT nr, const INT nc,REAL *a,INT *p);
/***************************************************************/
gcomplex *form_macro(input_grid *g){
  gcomplex *gc=malloc(sizeof(gcomplex));
  gc->ncells=1; //one macroelement only
  gc->n=g->dim;
  gc->nv=g->nv;
  gc->ne=g->ne;
  gc->x=g->x;
  gc->xe=g->xe;
  /* 
     the array entries of nodes[] are determined (not written yet) by
     an external function which constructs macroelements from
     edges(segments).
  */
  return gc;
}
/***********************************************************************/
INT *set_input_grid(input_grid *g)
{
  /* 
     reorders the vertices lexicographically, uses this ordering and
     the input to define the edges of the macroelement.  Every edge is
     put into a subset, i.e. two edges (i1,i2) and (j1,j2) are
     considered equivalent iff (i2-i1)=(j2-j1).  The number of
     divisions in an equivalent set of edges is taken to be the
     largest from the equivalence class.  OUTPUT array is a "dim"
     array and for each direction gives the number of partitions.
  */
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
  for (i=0;i<g->ne;i++){
    j=g->seg[3*i+1]-g->seg[3*i];
    g->seg[3*i+2]=p[j]; 
    //    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  for (i=0;i<g->ne;i++){
    if(g->seg[3*i]) continue;
    j=g->seg[3*i+1]-g->seg[3*i]-1;
    p[j]=g->seg[3*i+2]; 
    fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  p=realloc(p,g->dim*sizeof(INT)); // realloc to dimension g->dim
  for (i=0;i<g->dim;i++){
    fprintf(stdout,"\ndirection:%d; div=%d",i,p[i]);
  }
  return p;
}
/*********************************************************************/
void map2mac(scomplex *sc,cube2simp *c2s, input_grid *g)
{
  /* 
     maps a uniform grid from the n-dimensional cube to a hexagonal 
     macroelement given by its coordinates xmacro[1:nvcube*dim]
     xmac[nvcube][dim]
  */
  INT i,j,k1,k2,k1c,k2c,kf,dim=sc->n,dim1 =sc->n+1;
  INT nv=sc->nv,ns=sc->ns;
  INT ksys;
  REAL *xmac=g->x;  
  REAL *xhat = (REAL *)calloc(dim,sizeof(REAL));
  REAL *xemac=(REAL *)calloc(c2s->ne*dim,sizeof(REAL));
  // convert midpoints from polar to cartesian.
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->labels[k1]];
    k2c=g->systypes[g->labels[k2]];
    fprintf(stdout,"\nverts=(%d,%d); coord_sys=(%d,%d)",k1,k2,k1c,k2c);fflush(stdout);
    if(g->labels[k1]==g->labels[k2] && k1c==1){
      //use xhat as a temp array:
      xhat[0]=0.5*(xmac[k1*dim]+xmac[k2*dim]);// this is rho
      // take half angles;
      for(j=1;j<dim;j++) {
	xhat[j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
      }
      polar2cart(dim,xhat,xemac+(i*dim));
      // translate by adding the origin. 
      ksys=g->labels[k1];// k1c and k2c should be the same below. 
			      //      k2c=g->labels[k2];
      for(j=0;j<dim;j++) {
	xemac[i*dim+j]+=g->ox[ksys*dim+j];
      }
    }
  }
  // end of mid points in polar;
  // now convert all vertices in cartesian as well. 
  for(i=0;i<c2s->nvcube;i++){
    k1c=g->systypes[g->labels[i]];
    if(k1c==1){
      memcpy(xhat,xmac+i*dim,dim*sizeof(REAL));
      polar2cart(dim,xhat,xmac+(i*dim));
      //      translate
    }
    ksys=g->labels[i];
    for(j=0;j<dim;j++) {
      xmac[i*dim+j]+=g->ox[ksys*dim+j];
    }
  }
  // now everything is in cartesian, so midpoints that are
  // not yet attended to are just averages.
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->labels[k1]];
    k2c=g->systypes[g->labels[k2]];
    //skip all polar mid points
    if(g->labels[k1]==g->labels[k2] && k1c==1) continue;
    fprintf(stdout,"\ncart:verts=(%d,%d); coord_sys=(%d,%d)",k1,k2,k1c,k2c);fflush(stdout);
    for(j=0;j<dim;j++) {
      xemac[i*dim+j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
    }
  }
  r2c(c2s->nvcube,dim,sizeof(REAL),xmac); // we need xmac by rows here
  r2c(c2s->ne,dim,sizeof(REAL),xemac); // we need xemac (mid points of
				       // edges) also by rows
  for(kf=0;kf<sc->nv;kf++){
    for(i=0;i<dim;i++)xhat[i]=sc->x[kf*dim+i];
    for(i=0;i<dim;i++){
      //      fprintf(stdout,"coordinate:%d",i);
      sc->x[kf*dim+i]=interp8(c2s,xmac+i*c2s->nvcube,xemac+i*c2s->ne,xhat);
      //      sc->x[kf*dim+i]=interp4(c2s,xmac+i*c2s->nvcube,xhat);
    }
  }
  //  r2c(dim,c2s->nvcube,sizeof(REAL),xmac); // we need xmac by columns here
  //  r2c(dim,c2s->ne,sizeof(REAL),xemac); // we need xemac by rows agin
  //  if(xhat) free(xhat);
  //  if(xemac) free(xemac);
  return;
}
INT main(INT argc, char **argv)
{
  INT i=-1,j=-1,k=-1,kperm;
  input_grid *g=parse_input_grid("grid.input");
  INT dim=g->dim;
  cube2simp *c2s=cube2simplex(dim);
  INT ns=c2s->ns,dim1=c2s->n+1;
  /*------------------------------------------------------*/
  scomplex *sc;
  INT intype=0;
  /*------------------------------------------------------*/
  INT *nd=set_input_grid(g);  
  /*this can be used to generate grids in a different way, but not now:*/
  /* if(intype<-1){ */
  /*   for(i=0;i<dim;i++) ndd[i]=1; */
  /*   sc=umesh(dim,ndd,c2s,intype); */
  /*   unirefine(nd,sc);     */
  /* }else{ */
  /*   sc=umesh(dim,nd,c2s,intype); */
  /* } */
  /*GENERATE UNIFORM GRID: nodes in each direction: ND*/
  sc=umesh(dim,nd,c2s,intype);
  if(nd) free(nd);
  fprintf(stdout,"\nuniform mesh in dim=%d; vertices: %d, simplexes %d\n",dim,sc->nv,sc->ns);
  fprintf(stdout,"\nedges=%d",c2s->ne);
  INT k1,k2,j1,j2,l1,l2;
  gcomplex *macs=form_macro(g);
  map2mac(sc,c2s,g);
  input_grid_free(g); 
  fprintf(stdout,"\nDone.\n");
  //  haz_scomplex_print(sc,0,"HAHA");
  if(dim==2||dim==3) {
    vtkw("newmesh.vtu",sc,0,0,1.);
  }
  cube2simp_free(c2s);
  haz_scomplex_free(sc);
  return 0;
}
  
/***********************************************************/
