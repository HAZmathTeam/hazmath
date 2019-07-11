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
/********************************FINCTIONS:*********************/
void cube2simp_free(cube2simp *c2s);
INT reverse(void *arr,INT length, size_t elsize);
cube2simp *cube2simplex(INT dim);
scomplex *umesh(const INT dim, INT *nd, cube2simp *c2s, const INT intype);
void polar2cart(INT dim, REAL *px, REAL *cx);
REAL interp8(cube2simp *c2s, REAL *u, REAL *ue, REAL *xhat);
REAL interp4(cube2simp *c2s, REAL *u, REAL *xhat);
void unirefine(INT *nd,scomplex *sc);
/***********************************************************************/
INT *set_input_grid1(input_grid *g,cube2simp *c2s)
{
  /* 
     Every edge is put into a subset, i.e. two edges (v(i1),v(i2)) and
     (v(j1),v(j2)) are considered equivalent iff (i2-i1)=(j2-j1).  The
     number of divisions in an equivalent set of edges is taken to be
     the largest from the equivalence class.  OUTPUT array is a "dim"
     array and for each direction gives the number of partitions.
  */
  INT i,j,k,iri,ici,pmem;
  pmem=2*g->nv;
  if(pmem<2*g->ne) pmem=2*g->ne;
  if(pmem<2*g->nel) pmem=2*g->nel;
  if(pmem<2*g->nf) pmem=2*g->nf;
  INT *p=calloc(pmem,sizeof(INT));// permutation and inverse permutation;
  //
  for (i=0;i<g->ne;i++){
    iri=g->seg[3*i];
    ici=g->seg[3*i+1];
    if(iri<ici){
      g->seg[3*i]=iri;
      g->seg[3*i+1]=ici;
    } else {
      g->seg[3*i]=ici;
      g->seg[3*i+1]=iri;
    }
    /* set up divisions */
    j=g->seg[3*i+1]-g->seg[3*i]; // should be always positive;
    //    fprintf(stdout,"\n%%z123=%d:(%d,%d);%d",i,3*i,3*i+1,g0->seg[3*efound[i]+2]);
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
  /*ORDER*/
  ilexsort(g->ne, 3,g->seg,p);
  k=0;
  for (i=0;i<g->ne;i++){
    if(g->seg[3*i]) continue;
    j=g->seg[3*i+1]-g->seg[3*i]-1;
    p[k]=g->seg[3*i+2];
    k++;
    //  fprintf(stdout,"\n[%d,%d]:div=%d",g->seg[3*i],g->seg[3*i+1],g->seg[3*i+2]);
  }
  p=realloc(p,g->dim*sizeof(INT)); // realloc to dimension g->dim
  //  for (i=0;i<g->dim;i++){
  //    fprintf(stdout,"\ndirection:%d; div=%d",i,p[i]);
  //  }
  //  input_grid_print(g);
  //  print_full_mat_int(g->ne,3,g->seg,"med");
  //  print_full_mat_int(g->nf,(c2s->nvface+1),g->mfaces,"mf");
  //  print_full_mat_int(g->nel,(c2s->nvcube+1),g->mnodes,"mel");
  return p;
}
/***********************************************************************/
INT set_ndiv_edges(input_grid *g,		\
		   input_grid *g0,		\
		   cube2simp *c2s,		\
		   INT *efound,			\
		   INT *ffound,			\
		   INT **nd,
		   const INT iter)
{
  /* 
     For a given global input grid g0 creates local input grids for
     every macroelement and computes the divisions for it. It is used
     iteratively in macro_split to se the correct divisions for every
     macroelement. efound should be of length c2s->ne +g0->ne and
     ffound should be of length c2s->nf +g0->nf.  The input_grid g0
     should be all set, and the input_grid g should have all its
     scalar values set.
     nd is the array with the divisions, it must be g0->nel by c2s->n. 
  */
  INT kel0,i,j0,j1,swp,kel,ke,k0,k1,ndiv,pmem;
  INT nel0=g0->nel,nvcube=c2s->nvcube,nvface=c2s->nvface;
  INT *e0found=efound + c2s->ne;
  INT *f0found=ffound + c2s->nf;
  /*foe easier reference*/
  for(i=0;i<g0->ne;i++)e0found[i]=-1;
  for(i=0;i<g0->nf;i++)f0found[i]=-1;
  // make all divisions > 0
  for(ke=0;ke<g0->ne;ke++){
    ndiv=abs(g0->seg[3*ke+2]);
    if(ndiv==0)ndiv=1;
    g0->seg[3*ke+2]=ndiv;
  }
  for(ke=0;ke<g0->ne;ke++)      
    e0found[ke]=g0->seg[3*ke+2];
  for(kel0=0;kel0<nel0;kel0++){
    if((iter%2)) kel=nel0-kel0-1; else kel=kel0;
    // macroelement by macroelement try to find the edge divisions 
    for(i=0;i<c2s->ne;i++){
      g->seg[3*i]=c2s->edges[2*i];
      g->seg[3*i+1]=c2s->edges[2*i+1];
      g->seg[3*i+2]=-1;
      efound[i]=-1;
    }
    memcpy(g->mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
    for(i=0;i<c2s->ne;i++){
      j0=g->mnodes[c2s->edges[2*i]];
      j1=g->mnodes[c2s->edges[2*i+1]];
      if(j0>j1){swp=j0;j0=j1;j1=swp;}
      for(ke=0;ke<g0->ne;ke++){
	k0=g0->seg[3*ke];
	k1=g0->seg[3*ke+1];
	if((k0==j0)&&(k1==j1)){
	  g->seg[3*i+2]=g0->seg[3*ke+2];
	  efound[i]=ke;
	}
      }
      //	  fprintf(stdout,"\nElement:%d, edge=(%d,%d);",kel,j0,j1);
    }
      //    input_grid_print(g);
    nd[kel]=set_input_grid1(g,c2s);
    for(i=0;i<g->ne;i++){      
      if(efound[i]<0) continue;
      ke=efound[i];
      g0->seg[3*ke+2]=g->seg[3*i+2];
    }
    //    print_full_mat_int(g0->ne,3,g0->seg,"edglob");
  }
  INT chng=0;
  for(ke=0;ke<g0->ne;ke++){
    k1=abs(e0found[ke]-g0->seg[3*ke+2]);
    if(k1>chng)chng=k1;
  }
  //  fprintf(stderr,"\nchng=%d",chng);
  return chng;
}
/*******************************************************/
scomplex *macro_split(input_grid *g0,cube2simp *c2s)
{
  /* 
     From an input grid loops over the macroelements and sets up the
     divisions in every dimension. First makes the array of edges with
     their division consistent (the input can be inconsistent) grids
     each having a single macroelement.
  */
  input_grid *g;
  scomplex *sc;
  INT i,j0,j1,kel,ke,pmem;
  INT nel0=g0->nel,nvcube=c2s->nvcube,nvface=c2s->nvface;
  pmem=2*g0->nv;
  if(pmem<2*g0->ne) pmem=2*g0->ne;
  if(pmem<2*g0->nel) pmem=2*g0->nel;
  if(pmem<2*g0->nf) pmem=2*g0->nf;
  INT *p=calloc(pmem,sizeof(INT));
  ilexsort(g0->nel,(c2s->nvcube+1),g0->mnodes,p);
  ilexsort(g0->nf, (c2s->nvface+1),g0->mfaces,p);
  /*-------------------------------------------------------------------*/
  INT *efound=calloc(c2s->ne+g0->ne,sizeof(INT));
  INT *ffound=calloc(c2s->nf+g0->nf,sizeof(INT));
  /*-------------------------------------------------------------------*/
  INT **nd=calloc(g0->nel,sizeof(INT *));
  INT **nfcodes=calloc(g0->nel,sizeof(INT *));
  INT **neib=calloc(g0->nel,sizeof(INT *));
  for(i=0;i<g0->nel;i++){
    nd[i]=calloc(c2s->n,sizeof(INT)); /* to hold the number of
					 divisions in every coordinate
					 direction */
    nfcodes[i]=calloc(c2s->nf,sizeof(INT)); /* to hold the codes of
					       the faces */
    neib[i]=calloc(c2s->nf,sizeof(INT));/* to hold the element
					   neighboring list */
  }
  /*-------------------------------------------------------------------*/
  g=malloc(1*sizeof(input_grid));
  /**/
  g->title=g0->title;
  g->dgrid=g0->dgrid;
  g->fgrid=g0->fgrid;
  g->dvtu=g0->dvtu;
  g->fvtu=g0->fvtu;
  g->print_level=g0->print_level;
  g->ref_type=g0->ref_type;
  g->nref=g0->nref;
  g->err_stop=g0->err_stop;
  /**/
  g->dim=c2s->n;
  g->ncsys=g0->ncsys;
  g->nv=c2s->nvcube;
  g->nf=c2s->nf;
  g->ne=c2s->ne;
  g->nel=1;
  input_grid_arrays(g);
  /* reassign this as these are the same as g0 */
  free(g->systypes);   g->systypes=g0->systypes;
  free(g->syslabels);   g->syslabels=g0->syslabels;
  free(g->ox); g->ox=g0->ox;
  /* g->mnodes=(INT *)calloc(nvcube+1,sizeof(INT)); */
  /* g->mfaces=(INT *)calloc(nvface+1,sizeof(INT)); */
  /* g->seg=(INT *)calloc(3*g->ne,sizeof(INT)); */
  /* g->csysv=(INT *)calloc(g->nv,sizeof(INT)); */
  /* g->labelsv=(INT *)calloc(g->nv,sizeof(INT)); */
  /* g->bcodesv=(INT *)calloc(g->nv,sizeof(INT)); */
  /* g->xv=(REAL *)calloc(g->dim*g->nv,sizeof(REAL));  */
  /* g->xe=(REAL *)calloc(g->dim*g->ne,sizeof(REAL)); */
  INT chng=1,iter=0,maxiter=1024;  
  while(chng&&(iter<maxiter)){
    iter++;
    // make the divisions in g0->seg consistent;
    chng=set_ndiv_edges(g,g0,c2s,efound,ffound,nd,iter);
  }
  /* set the divisions on every edge now; since they are consistent we
   have: */
  if(set_ndiv_edges(g,g0,c2s,efound,ffound,nd,0)) {
    fprintf(stderr,"\n\n***ERR in %s: the divisions of the edges cannod be inconsistent during second call of set_ndiv_edges()\n\n",__FUNCTION__);
    exit(4);
  }
  // form macroelement neighboring list. Use transposition
  iCSRmat *el2v=malloc(1*sizeof(iCSRmat));
  el2v[0]=icsr_create(g0->nel,g0->nv,nvcube*g0->nel);
  //  iCSRmat el2v=icsr_create(g0->nel,g0->nv,nvcube*g0->nel);
  el2v->IA[0]=0;
  for(kel=0;kel<nel0;kel++){
    memcpy((el2v->JA+el2v->IA[kel]),(g0->mnodes+kel*(nvcube+1)),nvcube*sizeof(INT));
    el2v->IA[kel+1]=el2v->IA[kel]+nvcube;
  }
  for(i=0;i<el2v->IA[nel0];i++)
    el2v->val[i]=1;
  //  fprintf(stdout,"\nel2v=[");
  //  icsr_print_matlab_val(stdout,el2v);
  //  fprintf(stdout,"];");
  iCSRmat *v2el=malloc(1*sizeof(iCSRmat));
  icsr_trans(el2v,v2el);	
  //  fprintf(stdout,"\nv2el=[");
  //  icsr_print_matlab_val(stdout,v2el);
  //  fprintf(stdout,"];\n");
  iCSRmat *el2el=malloc(1*sizeof(iCSRmat));
  icsr_mxm(el2v,v2el,el2el);
  //  fprintf(stdout,"\nel2el0=[");
  //  icsr_print_matlab_val(stdout,el2el);
  //  fprintf(stdout,"];\n\n");
  icsr_free(el2v); icsr_free(v2el);
  // shrink the el2el matrix;
  INT jel;
  el2el->nnz=el2el->IA[0];
  for(kel=0;kel<nel0;kel++){    
    j0=el2el->IA[kel];
    j1=el2el->IA[kel+1];
    el2el->IA[kel]=el2el->nnz;
    for(ke=j0;ke<j1;ke++){
      jel=el2el->JA[ke];
      if((jel==kel) || (el2el->val[ke]!=nvface)) continue;
      el2el->JA[el2el->nnz]=el2el->JA[ke];
      el2el->val[el2el->nnz]=el2el->val[ke];
      el2el->nnz++;
    }
  }
  el2el->IA[nel0]=el2el->nnz;
  el2el->JA=realloc(el2el->JA,el2el->nnz*sizeof(INT));
  el2el->val=realloc(el2el->val,el2el->nnz*sizeof(INT));
  //    fprintf(stdout,"\nel2el=[");
  //    icsr_print_matlab_val(stdout,el2el);
  //    fprintf(stdout,"];\n\n");
  for(kel=0;kel<nel0;kel++){
    memcpy(g->mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
    for(i=0;i<g->nf;i++){
      for(ke=0;ke<nvface;ke++){
	g->mfaces[i*(nvface+1)+ke]=c2s->faces[i*nvface+ke];
      }
      g->mfaces[i*(nvface+1)+nvface]=-1;
    }
    //    print_full_mat_int(g->nf,(nvface+1),g->mfaces,"mf{1}"); fflush(stdout);
  }
  for(kel=0;kel<nel0;kel++){
    memcpy(g->mnodes,(g0->mnodes+kel*(nvcube+1)),(nvcube+1)*sizeof(INT));
    for(i=0;i<g->nv;i++){
      j0=g->mnodes[i];// vertex number (global)
      g->csysv[i]=g0->csysv[j0]; 
      g->labelsv[i]=g0->csysv[j0];
      memcpy((g->xv+i*g->dim),(g0->xv+j0*g0->dim),g->dim*sizeof(REAL));
    }
    //    print_full_mat(g->nv,g->dim,g->xv,"xv{1}"); fflush(stdout);
  }
  //////
  free(efound);
  free(ffound);
  for(i=0;i<g0->nel;i++) {
    free(nd[i]);
    free(nfcodes[i]);
    free(neib[i]);
  }
  free(nd);
  free(nfcodes);
  free(neib);
  exit(33);
  return sc;  
}
/*********************************************************************/
void map2mac(scomplex *sc,cube2simp *c2s, input_grid *g)
{
  /* 
     maps a uniform grid from the n-dimensional cube to a hexagonal 
     macroelement given by its coordinates xmacro[1:nvcube*dim]
     xmac[nvcube][dim]
  */
  INT i,j,k1,k2,k1c,k2c,kf,dim=sc->n;
  INT ksys;
  REAL *xmac=g->xv;  
  REAL *xhat = (REAL *)calloc(dim,sizeof(REAL));
  REAL *xemac=(REAL *)calloc(c2s->ne*dim,sizeof(REAL));
  // convert midpoints from polar to cartesian.
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->csysv[k1]];
    k2c=g->systypes[g->csysv[k2]];
    fprintf(stdout,"\nverts=(%d,%d); coord_sys=(%d,%d)",k1,k2,k1c,k2c);fflush(stdout);
    if(g->csysv[k1]==g->csysv[k2] && k1c==1){
      //use xhat as a temp array:
      xhat[0]=0.5*(xmac[k1*dim]+xmac[k2*dim]);// this is rho
      // take half angles;
      for(j=1;j<dim;j++) {
	xhat[j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
      }
      polar2cart(dim,xhat,xemac+(i*dim));
      // translate by adding the origin. 
      ksys=g->csysv[k1];// k1c and k2c should be the same below. 
			      //      k2c=g->csysv[k2];
      for(j=0;j<dim;j++) {
	xemac[i*dim+j]+=g->ox[ksys*dim+j];
      }
    }
  }
  // end of mid points in polar;
  // now convert all vertices in cartesian as well. 
  for(i=0;i<c2s->nvcube;i++){
    k1c=g->systypes[g->csysv[i]];
    if(k1c==1){
      memcpy(xhat,xmac+i*dim,dim*sizeof(REAL));
      polar2cart(dim,xhat,xmac+(i*dim));
      //      translate
    }
    ksys=g->csysv[i];
    for(j=0;j<dim;j++) {
      xmac[i*dim+j]+=g->ox[ksys*dim+j];
    }
  }
  // now everything is in cartesian, so midpoints that are
  // not yet attended to are just averages.
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->csysv[k1]];
    k2c=g->systypes[g->csysv[k2]];
    //skip all polar mid points
    if(g->csysv[k1]==g->csysv[k2] && k1c==1) continue;
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
  if(xhat) free(xhat);
  if(xemac) free(xemac);
  return;
}
INT main(INT argc, char **argv)
{
  //  INT j=-1;
  input_grid *g=parse_input_grid("grid.input");
  INT dim=g->dim;
  cube2simp *c2s=cube2simplex(dim);
  /*------------------------------------------------------*/
  scomplex *sc;
  INT intype=0;
  /*------------------------------------------------------*/
  macro_split(g,c2s);
  INT *nd=set_input_grid1(g,c2s);  
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
  fprintf(stdout,"\nGenerated a uniform mesh in dim=%d; vertices: %d, simplexes %d",dim,sc->nv,sc->ns);
  if(nd) free(nd);
  fprintf(stdout,"\nedges=%d",c2s->ne);
  //  INT k1,k2,j1,j2,l1,l2;
  //KEEP KEEP  gcomplex *macs=form_macro(g);
  fprintf(stdout,"\nMapping back to the macroelement...\n");
  map2mac(sc,c2s,g);
  input_grid_free(g); 
  fprintf(stdout,"\nDone.\n");
  //  haz_scomplex_print(sc,0,"HAHA");
  if(dim==2||dim==3) {
    fprintf(stdout,"Writing vtk file...\n");
    vtkw("newmesh.vtu",sc,0,0,1.);
  }
  cube2simp_free(c2s);
  haz_scomplex_free(sc);
  return 0;
}
/*********************EOF**********************************/
