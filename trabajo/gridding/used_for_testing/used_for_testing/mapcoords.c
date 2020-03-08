#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "hazmath.h"
#include "grid_defs.h"
#include "grid_params.h"
/*transforms [0,1]^d to a quad by bilinear transform what to do with
  the polar coords: implement later */
/************************************************************/
void zinterp(unigrid *ug,REAL *x, const INT nvert, INT *mask, const INT maskvalue)
/* here x is (ug->n+1) x nvert; the first ug->n columns are coords of
   indep. variables and the last coordinate is where the values are
   written. */
void trans0(const INT nv,const INT dim,		\
	    REAL *xmacro,REAL *x,INT *ib)
{
  /* 
     xmacro[] are coordinates of 2^d vertices of an n-dimensional
     parallelepiped. the array x[] is input and output. The mapping of
     vertices is as follows: if xo=min[x[j]] and xn=max[x[j]],
     j=0,(dim-1)then the parallelepiped defined by xo and xn is mapped
     to the parallelepiped defined by xmacro[]. x macro is ordered in
     increasing lexicographical order.
  */
{  
  INT i,j,mi1,k,l,kdim;
  REAL xoj,s;
  unsigned int bi,bi1;
  unsigned int *bits=ug->bits;
  INT nvcube=ug->nvcube;
  //    fprintf(stderr,"\nnvcube=%d\n",nvcube);
  INT *nd=ug->ndiv;
  INT *mo=(INT *)calloc(dim,sizeof(INT));
  REAL *u=(REAL *)calloc(nvcube,sizeof(REAL));
  REAL *xhat=(REAL *)calloc(dim,sizeof(REAL));
  /*NOTE: nd[] are the number of DIVISIONS (not number of vertices)
    for every direction */
  REAL eps0=1e-8;
  //  exit(12);
  for(j=0;j<dim;j++){
    dx[j]=(xn[j]-xo[j]);
  }
  
  for (i = 0; i<nvert;i++){
    for(j=0;j<dim;j++){ 
      xhat[j] = (x[dim*i+j]-xo[j])/dx[j];
      mo[j] = (INT )floor(xhat[j]-eps0);
      if(mo[j]<0 && fabs(xhat[j])<eps0*2.) {mo[j]=0;}
      else if(mo[j]<0) {mo[j]=0; continue;}
      if(mo[j]>nd[j]) {mo[j]=nd[j]; continue;}
      found=TRUE;
    }   
    if(!found){
      fprintf(stdout,"\nvertex=%i NOT found. x (xhat)=",i); fflush(stdout);
      for(j=0;j<dim;j++){ 
	fprintf(stdout," %12.6e (%10.5f)",x[dimbig*i+j],xhat[j]);
      }
      fprintf(stdout,"\n");
      //    }else{
      //      fprintf(stdout,"\nvertex=%i found.",i); fflush(stdout);
    }
    //    if(!found)continue;
    for(j=0;j<dim;j++){
      xoj=xo[j]+((REAL )mo[j])*dx[j];
      xhat[j] = (x[dim*i+j]-xoj)/dx[j];
    }
    /* for(j=0;j<dim;j++){ */
    /*   fprintf(stdout,"xo=%10.5e; dx=%12.3e ; xhat=%10.5f (%10.5f)\n ",xo[j],dx[j],xhat[j],x[dim*i+j]); fflush(stdout); */
    /* } */
    for kdim
    for(k = 0;k<nvcube;k++){
      kdim=k*dim;	
      for (j=0;j<dim; j++){mo[j]+=(INT )(!bits[kdim+j]);}
      /* return to the previous state */
      for (j=0;j<dim; j++){mo[j]-=(INT )(!bits[kdim+j]); }
      /* set the local values at the vertices */
      //      u[k]=ug->data[kf];
      u[k]=xmacro[dim*k+jj]
    }
    /* for(j=0;j<nvcube;j++){ */
    /*   fprintf(stdout,"u=%10.5e ",u[j]); fflush(stdout); */
    /* } */
    u=xmacro+
    x[i*dim+dim-1] = interp01(dim,bits,u,xhat); //already in hazmath
    /*    fprintf(stdout,"Err(%i)=%12.3e;\n",i+1,fi[i]-ff(dim,(x+dim*i)));*/
  }
  if(u) free(u);
  if(mo) free(mo);
  if(xhat) free(xhat);
  return;
}

