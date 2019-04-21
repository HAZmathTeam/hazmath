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
#ifndef DIM
#define DIM 2
#endif
/********************************FINCTIONS:*********************/
void cube2simp_free(cube2simp *c2s);
INT reverse(void *arr,INT length, size_t elsize);
cube2simp *cube2simplex(INT dim);
scomplex *umesh(const INT dim, INT *nd, cube2simp *c2s, const INT intype);
void polar2cart(INT dim, coordsystem *polar,REAL *px, REAL *cx);
//////////////////////////////////////////////////////////////////////
void unirefine(INT *nd,scomplex *sc);
/***************************************************************/
REAL interp4(cube2simp *c2s, REAL *u, REAL *xhat)
{
  /*INTerpolate d-linearly in d dimensions on the UNIT cube */
  INT dim=c2s->n,k,i,kdim;
  REAL phik;
  REAL s=0.;
  for(k = 0;k<c2s->nvcube;k++){
    kdim=k*dim;
    phik=1e0;
    for(i = 0;i<dim;++i){      
      if(!c2s->bits[kdim+i])
	phik*=(1e0-xhat[i]);
      else	
	phik*=xhat[i];
    }
    s+=u[k]*phik;
  }
  return s;
}
REAL interp8(cube2simp *c2s, REAL *u, REAL *ue, REAL *xhat)
{
  /*INTerpolate quadratically in d dimensions on the UNIT cube */
  INT dim=c2s->n,k,i,j,l,k1,k2,kdim;  
  REAL rl,phik,zmid,psimid;
  REAL s=0.;
  for(k = 0;k<c2s->nvcube;k++){
    kdim=k*dim;
    phik=1e0;
    for(i = 0;i<dim;++i){      
      if(!c2s->bits[kdim+i]){
	phik*=(1e0-xhat[i]);
	//	fprintf(stdout,"\nold:%d:vert:(1-x[%d])",k,i);
      } else {
	phik*=xhat[i];
	//	fprintf(stdout,"\nold:%d:vert:x[%d]",k,i);
      }
    }
    psimid=0.;zmid=0.;
    //    fprintf(stdout,"\nvertex=%d; normal=",k);
    for(i=0;i<dim;++i){
      if(!c2s->bits[kdim+i]){
	psimid-=xhat[i];
	zmid+=0.25e00;
	//	fprintf(stdout,"\n%d:vert:(-x[%d])",k,i);
      } else {
	psimid+=xhat[i];
	zmid-=0.75e00;
	//	fprintf(stdout,"\n%d:vert:x[%d]",k,i);
      }
    }
    //    fprintf(stdout,"\n%d:zmid=%5.1f",k,zmid);
    phik*=2e00*(zmid+psimid);
    s+=u[k]*phik;
  }
  REAL se,phie,xe1,xe2;
  se=0.;
  for(k=0;k<c2s->ne;k++){
    k2=c2s->edges[2*k];    
    k1=c2s->edges[2*k+1];
    fprintf(stdout,"\nmid=(%d,%d);",k1,k2);
    phie=1e0;
    for(i = 0;i<dim;++i){      
      if(!c2s->bits[k1*dim+i]){
	phie*=(1e0-xhat[i]);
	//	fprintf(stdout,"\n1:(1-x[%d])",i);
      } else {
	phie*=xhat[i];
	//	fprintf(stdout,"\n1:x[%d]",i);
      }
    }
    for(i = 0;i<dim;++i){      
      if(c2s->bits[k2*dim+i]==c2s->bits[k1*dim+i]) continue;
      if(!c2s->bits[k2*dim+i]){
	phie*=(1e0-xhat[i]);
	//	fprintf(stdout,"\n2:(1-x[%d])=%e",i,1.-xhat[i]);
      } else {
	phie*=xhat[i];
	//	fprintf(stdout,"\n2:x[%d]=%e",i,xhat[i]);
      }
    }
    phie*=4e00;
    fprintf(stdout,"\nedge=(%d,%d),coord=%e,phie=%e",k1,k2,ue[k],phie);
    se+=ue[k]*phie;
  }
  //  fprintf(stdout,"\n");
  fprintf(stdout,"\n***************** xhat=(%e,%e):%e",xhat[0],xhat[1],se);
  //  fprintf(stdout,"\ns=%e",s);
  fprintf(stdout,"\n");
  return (s+se);
}
void map2mac(scomplex *sc,cube2simp *c2s,REAL *xmac)
{
  /* 
     maps a uniform grid from the n-dimensional cube to a hexagonal 
     macroelement given by its coordinates xmacro[1:nvcube*dim]
     xmac[nvcube][dim]
  */
  INT i,j,k1,k2,kf,dim=sc->n,dim1 =sc->n+1;
  INT nv=sc->nv,ns=sc->ns;
  REAL *xhat = (REAL *)calloc(dim,sizeof(REAL));
  REAL *xemac=(REAL *)calloc(c2s->ne*dim,sizeof(REAL));
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    for(j=0;j<dim;j++) {
      if(k1<0){
	xemac[i*dim+j]=xmac[k1*dim+j]+0.25*xmac[k2*dim+j];
      } else {
	xemac[i*dim+j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
      }
    }
  }  
  /*********** PRINTING */
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    fprintf(stdout,"\nedge=(%d,%d); mid=(",k1,k2);
    for(j=0;j<dim-1;j++)fprintf(stdout,"%10.3e,",xemac[i*dim+j]);
    fprintf(stdout,"%10.3e); [",xemac[i*dim+dim-1]);
    for(j=0;j<dim-1;j++)fprintf(stdout,"%d,",c2s->bits[k1*dim+j]);
    fprintf(stdout,"%d]---[",c2s->bits[k1*dim+dim-1]);
    for(j=0;j<dim-1;j++) fprintf(stdout," %d ",c2s->bits[k2*dim+j]);
    fprintf(stdout,"%d]; (",c2s->bits[k2*dim+dim-1]);
    for(j=0;j<dim-1;j++)fprintf(stdout,"%5.1f,",xmac[k1*dim+j]);
    fprintf(stdout,"%5.1f)===(",xmac[k1*dim+dim-1]);
    for(j=0;j<dim-1;j++) fprintf(stdout," %5.1f ",xmac[k2*dim+j]);
    fprintf(stdout,"%5.1f)",xmac[k2*dim+dim-1]);
  }  
  fprintf(stdout,"\n\n");
  /*********** END PRINTING */
  r2c(c2s->nvcube,dim,sizeof(REAL),xmac); // we need xmac by rows here
  r2c(c2s->ne,dim,sizeof(REAL),xemac); // we need xemac (mid points of
				       // edges) also by rows
  for(kf=0;kf<sc->nv;kf++){
    for(i=0;i<dim;i++)xhat[i]=sc->x[kf*dim+i];
    for(i=0;i<dim;i++){
      fprintf(stdout,"coordinate:%d",i);
      sc->x[kf*dim+i]=interp8(c2s,xmac+i*c2s->nvcube,xemac+i*c2s->ne,xhat);
//      sc->x[kf*dim+i]=interp4(c2s,xmac+i*c2s->nvcube,xhat);
    }
  }
  r2c(dim,c2s->nvcube,sizeof(REAL),xmac); // we need xmac by columns here
  r2c(dim,c2s->ne,sizeof(REAL),xemac); // we need xemac by rows agin
  if(xhat) free(xhat);
  if(xemac) free(xemac);
  return;
}
INT main(INT argc, char **argv)
{
  INT i=-1,j=-1,k=-1,kperm,dim=DIM;
  cube2simp *c2s=cube2simplex(dim);
  INT ns=c2s->ns,dim1=c2s->n+1;
  INT *nd=(INT *)calloc(dim,sizeof(INT));  
  INT *ndd=(INT *)calloc(dim,sizeof(INT));  
  for(i=0;i<dim;i++) nd[i]=2 + 3*i;
  INT memx=dim;
  for(i=0;i<dim;i++) memx*=(nd[i]+1);
  /*------------------------------------------------------*/
  scomplex *sc;
  INT intype=0;
  /*------------------------------------------------------*/
  if(intype<-1){
    for(i=0;i<dim;i++) ndd[i]=1;
    sc=umesh(dim,ndd,c2s,intype);
    unirefine(nd,sc);    
  }else{
    sc=umesh(dim,nd,c2s,intype);
  }
  if(nd) free(nd);
  if(ndd) free(ndd);  
  fprintf(stdout,"\nuniform mesh in dim=%d; vertices: %d, simplexes %d\n",dim,sc->nv,sc->ns);
  fprintf(stdout,"\nedges=%d",c2s->ne);
  /* REAL xmac[24]={-1,-1.,-1.,				\ */
  /* 		 1.,-1,-1.,				\ */
  /* 		 -1., 1.,-1.,				\ */
  /* 		 1., 1.,-1.,				\ */
  /* 		 -1.,-1., 1.,				\ */
  /* 		 1.,-1., 1.,				\ */
  /* 		 -1., 1., 1.,				\ */
  /* 		 1., 1., 1.}; */
  
  INT k1,k2,j1,j2,l1,l2;
  /* if(dim==2){ */
  /*   xmac[0]=-1; xmac[1]=-1.; */
  /*   xmac[2]= 1.; xmac[3]=-1.; */
  /*   xmac[4]=-1.; xmac[5]= 1.; */
  /*   xmac[6]= 1.; xmac[7]= 1.;    */
  /*   for(i=8;i<24;i++){ */
  /*     xmac[i]=-1e20; */
  /*   }     */
  /* }   */
  input_grid *g=parse_input_grid("input.grid");
  //  
  input_grid_print(g);  
  map2mac(sc,c2s,g->x);
  input_grid_free(g);
  exit(129);
  //  haz_scomplex_print(sc,0,"HAHA");
  if(dim==2||dim==3) {
    vtkw("newmesh.vtu",sc,0,0,1.);
  }
  cube2simp_free(c2s);
  haz_scomplex_free(sc);
  return 0;
}
  
