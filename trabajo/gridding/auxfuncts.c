#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "hazmath.h"
#include "grid_defs.h"
/* 
   file with auxiliary functions mainly supporting interpolation in
   R^n (from uniform grid to scattered data); general not necessarily
   related to hydro app. Even if the include
   above is commented out, this should compile.
*/
#ifndef REAL
#define REAL double
#endif
#ifndef INT
#define INT int
#endif
#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE  0
#endif
#ifndef nula
#define nula  -1
#endif

/* INT chkn(INT n, const INT nmin, const INT nmax) */
/* { */
/*   INT nnew=n; */
/*   /\*  */
/*      checks if n is in the closed interval [nmin, nmax] and if */
/*      not, changes the value of n to be the closest end of the */
/*      interval  */
/*   *\/ */
/*   if(n > nmax) { */
/*     fprintf(stdout, "\n Iput value too large: %d ; Changing to the max allowed: %d\n",n,nmax); */
/*     nnew=nmax; */
/*   } else if(n < nmin) { */
/*     fprintf(stdout, "\ninput value too small: %d ; Changing to the min allowed: %d\n",n,nmin); */
/*     nnew=nmin; */
/*   } */
/*   return nnew; */
/* } */

/********************************************************************/
/* interpolates using data from lattice grid; the coordinates of the
   points to interpolate at are given by x which is assumed to be of
   dimension ug->n+1; the interpolated values are stored as last column
   in x (which should be available for this). Anything outside the
   lattice is moved to the closest lattice boundary in L_1 and
   interpolated too.  */
void zinterp(unigrid *ug,REAL *x, const INT nvert, INT *mask, const INT maskvalue)
/* here x is (ug->n+1) x nvert; the first ug->n columns are coords of
   indep. variables and the last coordinate is where the values are
   written. */
{  
  INT dim = ug->n,dimbig=ug->n+1;
  INT i,j,mi1,k,kdim,kf;
  REAL xoj,s;
  unsigned int bi,bi1;
  REAL *xo = ug->xo, *xn = ug->xn, *dx = ug->dx;
  unsigned int *bits=ug->bits;
  INT nvcube=ug->nvcube;
  //    fprintf(stderr,"\nnvcube=%d\n",nvcube);
  INT *nd=ug->ndiv;
  INT *mo=(INT *)calloc(dim,sizeof(INT));
  REAL *u=(REAL *)calloc(nvcube,sizeof(REAL));
  REAL *xhat=(REAL *)calloc(dim,sizeof(REAL));
  /*NOTE: nd[] are the number of DIVISIONS (not number of vertices)
    for every direction */
  unsigned int found;
  REAL eps0=1e-8;
  //  exit(12);
  for (i = 0; i<nvert;i++){
    if(mask != NULL) {
      if(mask[i]!=maskvalue) continue; /* do not mask if not on maskvalue*/
    }
    found=FALSE;
    for(j=0;j<dim;j++){ 
      xhat[j] = (x[dimbig*i+j]-xo[j])/dx[j];
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
      xhat[j] = (x[dimbig*i+j]-xoj)/dx[j];
    }
    /* for(j=0;j<dim;j++){ */
    /*   fprintf(stdout,"xo=%10.5e; dx=%12.3e ; xhat=%10.5f (%10.5f)\n ",xo[j],dx[j],xhat[j],x[dim*i+j]); fflush(stdout); */
    /* } */
    for(k = 0;k<nvcube;k++){
      kdim=k*dim;	
      for (j=0;j<dim; j++){mo[j]+=(INT )(!bits[kdim+j]);}
      /* compute the global index */
      kf=num_lattice(mo,dim,nd);
      /* return to the previous state */
      for (j=0;j<dim; j++){mo[j]-=(INT )(!bits[kdim+j]); }
      /* set the local values at the vertices */
      u[k]=ug->data[kf];
    }
    /* for(j=0;j<nvcube;j++){ */
    /*   fprintf(stdout,"u=%10.5e ",u[j]); fflush(stdout); */
    /* } */
    x[i*dimbig+dimbig-1] = interp01(dim,bits,u,xhat); //already in hazmath
    /*    fprintf(stdout,"Err(%i)=%12.3e;\n",i+1,fi[i]-ff(dim,(x+dim*i)));*/
  }
  if(u) free(u);
  if(mo) free(mo);
  if(xhat) free(xhat);
  return;
}
void swne(INT n, INT nv, REAL *xo, REAL *xn, REAL *x)
{
  /*  finds sw nd ne corners of the rectangle containing x*/
  INT i,j;
  for(i = 0;i<n;i++){
    xo[i]=x[i];
    xn[i]=xo[i];
    for (j = 1;j<nv;j++){
      if(x[j*n+i]<xo[i]) xo[i]=x[j*n+i];
      if(x[j*n+i]>xn[i]) xn[i]=x[j*n+i];
    }
  }
  return;
}
void scalec(const INT flag, INT n, INT nv, REAL *xo, REAL *x, REAL *scale)
{
  /*maps x as described below; overwrites the input*/ 
  INT i,j;
  if (flag){
    /* FRW shift and then scale x<--(x-xo)*scale */
    for(i = 0;i<n-1;i++){
      for (j = 0;j<nv;j++)
	x[j*n+i]=(x[j*n+i]-xo[i])*scale[i];
    }
    /* treat z as the axis points down so it only scales */
    i=n-1;
    for (j = 0;j<nv;j++){
      x[j*n+i]=(x[j*n+i])*scale[i];
    }
  } else {
    /* BCK: unscale and then shift x<--x/scale+xo */
    REAL s;
    for(i = 0;i<n;i++){
      s=1./scale[i];
      for (j = 0;j<nv;j++)
	x[j*n+i]=x[j*n+i]*s+xo[i];
    }
  }
  return;
}
void scaleg(INT n, INT nv, REAL *xo, REAL *xn, REAL *x, REAL *scale)
{
  /*geo scale: maps x as described below with z (or y in 2d)-axis
    pointing down; overwrites the input*/ 
  INT i,j;
  for(i = 0;i<n-1;i++){
    for (j = 0;j<nv;j++)
  /* FRW ;shift and then scale x<--(x-xo)*scale */
      //      x[j*n+i]=(x[j*n+i]-xo[i])*scale[i];
  /* ATTN: no shift; keep origin the same x<--xo+(x-xo)*scale */
      x[j*n+i]=xo[i]+(x[j*n+i]-xo[i])*scale[i];
  }
  /* treat z as the axis points down so it only scales */
  i=n-1;
  for (j = 0;j<nv;j++){
    x[j*n+i]=(x[j*n+i]-xn[i])*scale[i];
    //    x[j*n+i]=(x[j*n+i]-xn[i])*scale[i]*95310./2000.;
  }
  return;
}
void mapcoords(REAL *xp, INT np, INT nbig,			\
	       INT n, REAL *xo0,REAL *xn0,REAL *xo1,REAL *xn1)
{
  REAL a,b,d;
  INT j,k;
  /* 
    
    maps the coordinates xp using the map from a rectangle xo0,xn0 to
    the rectangle defined by xo1,xn1 in Rn.  xp is (np by nbig)
    stored by rows and the first n columns are overwritten by the
    result. In another words, with nbig = 3 and n = 2 it transforms
    (xp,yp,zp) to (xpnew,ypnew,zp).
*/  
  for(j = 0;j<n;j++){
    d=1e0/(xn0[j]-xo0[j]);
    a=(xn1[j]-xo1[j])*d;
    b=(xo1[j]*xn0[j]-xn1[j]*xo0[j])*d;
    for(k = 0;k<np;k++) {
      xp[k*nbig+j]=a*xp[k*nbig+j]+b;
    }
  }
  /* fprintf(stdout,"\nMapping %d coords in %s:",np,__FUNCTION__); */
  /* print_full_mat(nbig,1,xo0,"xo0"); */
  /* print_full_mat(nbig,1,xn0,"xn0"); */
  /* print_full_mat(nbig,1,xo1,"xo1"); */
  /* print_full_mat(nbig,1,xn1,"xn1"); */
  return;
}
/***********************************************/
void vert_layer(INT nlayers,scomplex *sc,INT *mask, INT *nzptl)
{
  INT maxflag,node,dj,i,j,k,iflag,i1;
  INT dim=sc->n,dim1=dim+1,mi=-10;
  /*
    match vertices to layers: we have to asign to every vertex as
    well as to every simplex the layer number. For vertices this
    should be a flag from 1 to nlayers+1; for simplices this is
    from 1 to nlayers. on output from the meshing every simplex
    has a flag which equals the "nz" coordinate of its supporting
    plane.
  */
  INT *iwrk = (INT *)calloc(sc->nv,sizeof(INT));
  maxflag=-1;
  for(j=0;j<sc->ns;j++){
    i1=sc->flags[j]+1;
    if(maxflag<i1) maxflag=i1;
  }
  for(i=0;i<sc->nv;i++){mask[i]=-1;}
  for(i=0;i<sc->nv;i++){iwrk[i]=maxflag+2;}
  for(j=0;j<sc->ns;j++){
    dj=dim1*j;       
    iflag=sc->flags[j];
    i1=iflag+1;
    for(i=0;i<dim1;i++){
      node=sc->nodes[dj+i];
      //      if(mask[node]<(iflag+1)) {mask[node]=iflag+1;}
      if(mask[node]<(i1)) {mask[node]=i1;}
      if(iwrk[node]>(i1)) {iwrk[node]=i1;}
    }
  }
  for(i=0;i<sc->nv;i++){
    mask[i]=(mask[i]+iwrk[i]);
    if(mask[i]/2>=maxflag) mask[i]/=2;
    else mask[i]=(mask[i]-1)/2;
  }
  if(iwrk)free(iwrk);
  //
  for(i=0;i<sc->nv;i++){
    mask[i] = maxflag-mask[i]+1;
  }
  for(i=0;i<sc->ns;i++){ 
    sc->flags[i] = maxflag-sc->flags[i];
  }
  for(i=0;i<sc->ns;i++){
    mi=sc->flags[i];
    for(j=0;j<nlayers;j++){
      if(mi<=nzptl[j+1] && mi>nzptl[j]){
  	sc->flags[i]=j;
      }
    }
  }
  /* for(i=0;i<sc->nv;i++){  */
  /*   fprintf(stdout,"\n%d,z=%f,mask=%d",i,sc->x[i*dim+dim-1],mask[i]);  */
  /* }  */
  /*****************************************************/
  /*THIS BELOW IS NOT NEEDED */
  /* fprintf(stdout,"\nVERTICES\n"); */
  /* for(i=0;i<sc->nv;i++){ */
  /*   fprintf(stdout,"\nvert=%d; mask=%d, nzptl=",i,mask[i]); */
  /*   mi=mask[i]; */
  /*   for(j=1;j<nlayers+1;j++){ */
  /*     fprintf(stdout,"(%d:%d) ",nzptl[j-1],nzptl[j]); */
  /*     if(mi<=nzptl[j]+1 && mi>nzptl[j-1]+1){ */
  /* 	mask[i]=j; */
  /*     } */
  /*   } */
  /*   fprintf(stdout,"\nafter: el=%d; mask=%d",i,mask[i]); */
  /* } */
  /* fprintf(stdout,"\nELEMENTS\n"); */
  /* for(i=0;i<sc->ns;i++){ */
  /*   fprintf(stdout,"\nel=%d; flag=%d, nzptl=",i,sc->flags[i]); */
  /*   for(j=0;j<nlayers;j++){ */
  /*     fprintf(stdout,"(%d:%d) ",nzptl[j],nzptl[j+1]); */
  /*     } */
  /*   } */
  /*   fprintf(stdout,"\nafter: el=%d; flag=%d",i,sc->flags[i]); */
  /* } */
  return;
}


