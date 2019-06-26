/*! \file src/amr/unigrid.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *   \note routines to initialize construct and destruct uniform grids in
 *  d-dimensions for general d.
 *  \note interpolates data on uniform grids
*/
#include "hazmath.h"
void coord_lattice(INT *m,const INT dim,			\
		   const INT kf, const INT nall, const INT *nd)
{
  INT i,j,k;
  /*
    given a global number kf on a lattice grid with lexicographical
    ordering, this returns the n-tuple of latice coordinates,
    m[i],i=0:n-1). nall is the total number of vertices in the
    lattice, nd[i] is the number of divisions in every direction;
  */
  j=nall;  k=kf;
  for(i=dim;i>0;i--){
    j=j/(nd[i-1]+1);  m[i-1] = k/j;   k=k-m[i-1]*j;
  }
  return;
}
INT num_lattice(INT *m,const INT dim,INT *nd)
{
  INT i,j,k,kf;
  /*
    given the lattice coordinates m[i],i=0:n-1 of a vertex, in
    dimension n, finds the global number kf of a vertex on a
    lexicographically ordered lattice grid. divisions in each
    directions are stored in nd[].
  */
  kf=m[dim-1];
  //      kf = m[dim-1];
  for (i=dim-1; i>0; i--){
    kf=kf*(nd[i-1]+1)+m[i-1];
  }
  return kf;
}

void binary1(const INT dim, unsigned int *bits, INT *nvloc)
{
  // coordinates of the vertices of the unit dim-cube as arrays of 0/1 
  INT i,k,kdim=-10,nbits=dim-1;
  *nvloc = (1 << dim);
  for(k = 0;k<(*nvloc);k++){
    kdim=k*dim;
    for (i=nbits ; i >=0; --i){
      bits[kdim+i] = (unsigned int )(!(k >> i & 1));
    }
  }
  return;
}

unigrid *ugrid_init(INT n, INT *nd, REAL *xo, REAL *xn)
{
  INT j,i,k;
  unigrid *ug=(unigrid *)malloc(sizeof(unigrid));
  /* fprintf(stdout,"\n"); */
  /* for(j=0;j<n;j++){ */
  /*   fprintf(stdout,"%i ",nd[j]); */
  /* } */
  /* fprintf(stdout,"\n"); */
  /* for(j=0;j<n;j++){ */
  /*   fprintf(stdout,"%f %f\n",xo[j],xn[j]); */
  /* } */
  ug->n=n; 
  /* fprintf(stdout,"n=%i\n",ug->n);fflush(stdout); */
  ug->nvcube=(1<<n); 
  /* fprintf(stdout,"nvcube=%i\n",ug->nvcube);fflush(stdout); */
  ug->ndiv=nd; 
  ug->bits=(unsigned int *)calloc(n*(ug->nvcube),sizeof(unsigned int));
  binary1(n,ug->bits,&(ug->nvcube));
  /* fprintf(stdout,"\n"); */
  /* for(k = 0;k<ug->nvcube;k++){ */
  /*   for(i = 0;i<n;++i){ */
  /*     fprintf(stdout,"%i ",ug->bits[k*n+i]); fflush(stdout); */
  /*   } */
  /*   fprintf(stdout,"\n"); */
  /* } */
 /* */
 if(!xo || !xn){
   ug->xo=(REAL *)calloc(n,sizeof(REAL)); /* coordinates of the origin
					     xo[dim]*/
   ug->xn=(REAL *)calloc(n,sizeof(REAL)); /* coordinates of the max
					     corner(NE in 2D)
					     xn[dim]*/   
   for(j=0;j<n;j++){
     ug->xo[j]=0.;
     ug->xn[j]=1.;
   }
 } else {
   ug->xo=xo;
   ug->xn=xn;
 }
 ug->dx = (REAL *)calloc(n,sizeof(REAL));
 ug->nall=1;
 for(j=0;j<n;j++){
   ug->dx[j]=(ug->xn[j]-ug->xo[j])/((REAL )ug->ndiv[j]);
   ug->nall *= (ug->ndiv[j]+1);
   /* fprintf(stdout,"\ndx[%i]=%f\n",j,ug->dx[j]); */
   /* fprintf(stdout,"\n(%i)=%i: \n",(ug->ndiv[j] + 1),ug->nall); */
 }
 /* allocate space for the data */
 ug->data=(REAL *)calloc(ug->nall,sizeof(REAL));
 if(!ug->data) {
   fprintf(stdout,"\n***ERROR: could not allocate memory for a uniform grid with data in %s\n",__FUNCTION__);fflush(stdout);
   exit(7);
 }
 ug->ugtype = 0; 
 return ug; 
}
void ugrid_free(unigrid *ug)
{
  /* fprintf(stdout,"\n%li %li %li %li %li\n", */
  /* 	  (ug->ndiv), \ */
  /* 	  (ug->bits), \ */
  /* 	  (ug->xo),   \ */
  /* 	  (ug->xn),   \ */
  /* 	  (ug->dx));fflush(stdout); */
    //
  //  if(ug->ndiv) free(ug->ndiv);
  if(ug->bits) free(ug->bits);
  // if(ug->xo) free(ug->xo);
  // if(ug->xn) free(ug->xn);
  if(ug->dx) free(ug->dx);
  if(ug->data) free(ug->data);
  ug->data=NULL;
  if(ug) free(ug);
  ug=NULL;
  return; 
}
void ugrid_transform(const INT n,unigrid *ug,	\
		     REAL *xodst, REAL *xndst)
{
  /* transform the unform grid ug to be over a parallelepiped with
     corners xodst and xndst */
  INT j;
  //  REAL a,b,dxsrc,dxdst;
  ug->xo=xodst;
  ug->xn=xndst;
  for(j=0;j<n;j++){
    /*    
	  dxsrc=ug->xn[j]-ug->xo[j];
	  dxdst=xndst[j]-xodst[j];
	  a=dxdst/dxsrc; b=(xodst[j]*ug->xn[j]-xndst[j]*ug->xo[j])/dxsrc;
    */
    ug->dx[j]=(ug->xn[j]-ug->xo[j])/((REAL )ug->ndiv[j]);
  }
  return;
}
/* general transform: maps data */
void data_transform(const INT nv, const int m,			\
		     REAL *data, REAL *xodst, REAL *xndst)
{
  /* transform each column of a REAL data given in (nv x m) matrix to
     the intervals xodst[j],xndst[j], j=1:m*/
  INT i,j;
  REAL *xo=(REAL *)calloc(2*m,sizeof(REAL));
  REAL *xn=xo+m;
  REAL dxsrc,a,b;
  for(j=0;j<m;j++){
    xn[j]=xo[j]=data[j];
    for(i=1;i<nv;i++){
      if(data[m*i+j]<xo[j]) xo[j]=data[m*i+j];
      if(data[m*i+j]>xn[j]) xn[j]=data[m*i+j];
    }
  }  
  for(j=0;j<m;j++){
    dxsrc=xn[j]-xo[j];
    a=(xndst[j]-xodst[j])/dxsrc;
    b=(xodst[j]*xn[j]-xndst[j]*xo[j])/dxsrc;    
    for(i=0;i<nv;i++){
      data[m*i+j] = a*data[m*i+j] +b;
    }
  }
  if(xo) free(xo);
  return; 
}
static REAL ff(const INT dim, REAL *x){
  /* function to interpolate (only for testing) */
  switch(dim){
  case 2:
    return 20e0+ 5.*x[0]-4.*x[1] + 17.*x[0]*x[1];
  case 3:
    return 20e0+5.*x[0]-4.*x[1]-0.3456*x[2]+x[0]*x[1];
  case 4:
    return 20e0+5.*x[0]-4.*x[1]-0.3456*x[2]-x[3];
  default:
    return 1e0;
  }
}
static void fdata(unigrid *ug)
{
  INT dim=ug->n,k,i,j,kf,nall=ug->nall;
  INT *nd=ug->ndiv;
  REAL *dx = ug->dx, *xo=ug->xo, *xn=ug->xn;
  REAL *x = (REAL *)calloc(dim,sizeof(REAL));
  INT *m = (INT *)calloc(dim,sizeof(INT));
  for (kf=0;kf<nall;kf++){
    //    j=floor(((double )nall)/((double )nd[dim-1]));
    coord_lattice(m,dim,kf,ug->nall,ug->ndiv);
    for(i=0;i<dim;i++){
      x[i]=xo[i]+m[i]*dx[i];
    }
    ug->data[kf]=ff(dim,x);
  }
  if(m) free(m);
  if(x) free(x);
  return;
}
REAL interp01(const INT dim,unsigned int *bits,REAL *u,	\
	      REAL *xhat)
{
  /*INTerpolate d-linearly in d dimensions on the UNIT cube */
  INT k,i,kdim,nvloc=(1<<dim);
  REAL phik;
  REAL s=0.;
  for(k = 0;k<nvloc;k++){
    kdim=k*dim;
    phik=1e0;
    for(i = 0;i<dim;++i){      
      if(bits[kdim+i])
	phik*=(1e0-xhat[i]);
      else	
	phik*=xhat[i];
    }
    s+=u[k]*phik;
  }
  return s;
}
/********************************************************************/
/* interpolates using data from lattice grid. Anything outside the
   lattice is moved to the closest lattice boundary in L_1 and
   interpolated too.  */
void interp1(const INT dimbig, REAL *fi, unigrid *ug,	\
	     REAL *x, const INT nvert, INT *mask)
/* here x is (dimbig) x nvert; but only first dim columns of it are used */
{  
  /*  find the values of ug->data at x by interpolation.  ug->data is
  given on a uniform grid.
  */
  INT dim = ug->n,i,j,mi1,k,kdim,kf;
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
  for (i = 0; i<nvert;i++){
    if(mask != NULL) {
      if(mask[i]) continue; /* mask false means do not mask */
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
    fi[i] = interp01(dim,bits,u,xhat);
    /*    fprintf(stdout,"Err(%i)=%12.3e;\n",i+1,fi[i]-ff(dim,(x+dim*i)));*/
  }
  if(u) free(u);
  if(mo) free(mo);
  if(xhat) free(xhat);
  return;
}
REAL interp02(const INT dim, uint *bits, REAL *u,	\
	      REAL *xhat)
{
  /*Interpolate d-linearly in d dimensions on the UNIT cube */
  // given a point with coordinates xhat[j],j=0:dim-1 in the unit
  // cube. and function values u[k] at the vertices of the unit cube
  // (k=0:2^{dim}-1)  this returns the value of the interpolant at xhat[]. 
  INT k,i,kdim,nvloc=(1<<dim);
  REAL phik;
  REAL s=0.;
  for(k = 0;k<nvloc;k++){
    kdim=k*dim;
    phik=1e0;
    for(i = 0;i<dim;++i){      
      if(bits[kdim+i])
	phik*=(1e0-xhat[i]);
      else	
	phik*=xhat[i];
    }
    s+=u[k]*phik;
  }
  return s;
}


/********************************************************************/
/* interpolates using data from lattice grid. Anything outside the
   lattice is moved to the closest lattice boundary in L_1 and
   interpolated too.  */
void interp2(REAL *fi, unigrid ug, REAL *x, const INT nvert, INT *mask)
{  
  /*  find the values of ug->data at array x by interpolation.  ug->data is
      given on a uniform grid.
      fi is vector to store the interpolation results 
  */
  INT dim = ug.n,i,j,mi1,k,kdim,kf;
  REAL xoj,s;
  uint bi,bi1;
  REAL *xo = ug.xo, *xn = ug.xn, *dx = ug.dx;
  uint *bits=ug.bits;
  INT nvcube=ug.nvcube;
  //    fprintf(stderr,"\nnvcube=%d\n",nvcube);
  INT *nd=ug.ndiv;
  INT *mo=(INT *)calloc(dim,sizeof(INT));
  REAL *u=(REAL *)calloc(nvcube,sizeof(REAL));
  REAL *xhat=(REAL *)calloc(dim,sizeof(REAL));
  /*NOTE: nd[] are the number of DIVISIONS (not number of vertices)
    for every direction */
  uint found;
  REAL eps0=1e-8;
  for (i = 0; i<nvert;i++){
    if(mask != NULL) {
      if(mask[i]) continue; /* mask false means do not mask */
    }
    found=FALSE;
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
	fprintf(stdout," %12.6e (%10.5f)",x[dim*i+j],xhat[j]);
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
    for(k = 0;k<nvcube;k++){
      kdim=k*dim;	
      for (j=0;j<dim; j++){mo[j]+=(INT )(!bits[kdim+j]);}
      /* compute the global index */
      kf=num_lattice(mo,dim,nd);
      /* return to the previous state */
      for (j=0;j<dim; j++){mo[j]-=(INT )(!bits[kdim+j]); }
      /* set the local values at the vertices */
      u[k]=ug.data[kf];
    }
    /* for(j=0;j<nvcube;j++){ */
    /*   fprintf(stdout,"u=%10.5e ",u[j]); fflush(stdout); */
    /* } */
    fi[i] = interp01(dim,bits,u,xhat);
    /*    fprintf(stdout,"Err(%i)=%12.3e;\n",i+1,fi[i]-ff(dim,(x+dim*i)));*/
  }
  if(u) free(u);
  if(mo) free(mo);
  if(xhat) free(xhat);
  return;
}
static REAL interp0(const INT dim,unsigned int *bits,REAL *u,	\
	    REAL *xo, REAL *dx, REAL *xstar)
{
  /*interpolate d-linearly in d dimensions on a single parallelepiped with
    one vertex given by xo[.] and edge lengths = dx[k], k=1:d*/
  INT k,i,kdim,nvloc=(1<<dim);
  REAL phik;
  REAL s=0.;
  for(k = 0;k<nvloc;k++){
    kdim=k*dim;
    phik=1e0;
    for(i = 0;i<dim;++i){      
      if(bits[kdim+i])
	phik*=(1e0-(xstar[i]-xo[i])/dx[i]);
      else	
	phik*=((xstar[i]-xo[i])/dx[i]);
    }
    s+=u[k]*phik;
  }
  return s;
}
static void fdataf(unigrid *ug)
{
  /*on a uniform grid ug fills a ug->data[] array with values computer with a function ff */ 
  INT dim=ug->n,k,i,j,kf,nall=ug->nall;
  INT *nd=ug->ndiv;
  REAL *dx = ug->dx, *xo=ug->xo, *xn=ug->xn;
  REAL *x = (REAL *)calloc(dim,sizeof(REAL));
  INT *m = (INT *)calloc(dim,sizeof(INT));
  for (kf=0;kf<nall;kf++){
    //    j=floor(((double )nall)/((double )nd[dim-1]));
    coord_lattice(m,dim,kf,ug->nall,ug->ndiv);
    for(i=0;i<dim;i++){
      x[i]=xo[i]+m[i]*dx[i];
    }
    ug->data[kf]=ff(dim,x);
  }
  if(m) free(m);
  if(x) free(x);
  return;
}
/************************ example of main.c to interpolate */
/* INT main(INT argc, char *argv[]) */
/* { */
/*   INT j,i; */
/*   INT dim = 3; */

/*   /\*TESTING*\/ */
/*   //  REAL xo[4]={0.,-1.,1.,1.};/\*SW corner coordinates*\/ */
/*   //  REAL xn[4]={1.,2.,3.,5.}; /\*NE corner coordinates*\/ */

/*   REAL xo[4]={0.,-1.,1.};/\*SW corner coordinates*\/ */
/*   REAL xn[4]={1.,2.,3.}; /\*NE corner coordinates*\/ */

/*   INT nd[4]={3,4,5}; /\*number of divisions in every direction*\/ */
/*   INT nvert=3; /\*number of points to interpolate. coordinates below *\/  */
/*   REAL x[9]={0.123,0.,1.1,				\ */
/*   	     0.51,1.2345,1.5,				\ */
/*   	     0.1,1.8,2.5}; */
  
/*   unigrid *ug=ugrid_init(dim,nd, xo, xn); */

/*   fdataf(ug);  /\*calculate ff on uniform grid*\/ */
  
/*   REAL *fi=(REAL *)calloc(nvert, sizeof(REAL));  */
/*   interp1(fi,*ug,x,nvert,NULL); */
  
/*   // compute max norm error */
/*   REAL err0=-1., diff=-10.; */
/*   for(j=0;j<nvert;j++){ */
/*     diff=fabs(ff(dim,(x + j*dim))-fi[j]); */
/*     if(err0 < diff) err0=diff;  //err0 stores the largest difference */
/*   } */
/*   fprintf(stdout,"ERR=%e\n",err0); */

/*   /\*END TESTING*\/ */
/*   if(fi) free(fi); */
/*   //if(ug)  ugrid_free(ug); //THIS CAUSES ERROR */
/*   // if(ug) fprintf(stdout,"ug is here"); */
  
/*   return 0; */
/* } */
