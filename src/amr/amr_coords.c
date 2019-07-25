/*! \file src/amr/amr_coords.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note containing all routines for polar to cartesian; cartesian to
 *  polar, etc. It also has the routing that maps "quadratically" a macroelement. 
 *  Works in dimension n. 
 *
 * \todo add cyllindrical and other coords. (Ludmil)

 */
#include "hazmath.h"
/**********************************************************************/
static long double pi=4e00*atanl(1e00);
/**********************************************************************/
REAL deg2rad(REAL alpha_deg)
{
  // degrees to radians;
  return (REAL ) ((long double )alpha_deg)*(pi/180e00);
}
/**********************************************************************/
REAL zero_twopi(REAL alpha)
{
  // map alpha in radians to [0,2*pi]
  long double angle=(long double )alpha;
  angle=atan2l(sinl(angle),cosl(angle));
  if(fabsl(angle)<1e-15) angle=0e0;
  if(angle<0e0)
    return (REAL )(angle+pi+pi);
  return (REAL )angle;
}
REAL zero_twopi_deg(REAL alpha_deg)
{
  // map alpha in degrees to [0,2pi)
  return zero_twopi(deg2rad(alpha_deg));
}
/**********************************************************************/
static void coord_perm(SHORT type, INT n,void *x, size_t elsize)
{
  /*
    permutes coordinate arrays (or any array whose elements are of size): 
    if type=1: x[0],x[1]...x[n-1]--> x[1]...x[n-1],x[0] 
    if type=0(inverse): y[0],y[1]...y[n-1]--> y[n-1],y[0]...y[n-2] 
  */
  INT i;  
  void *xx0n=(void *)calloc(1,elsize*sizeof(void));
  if(type){
    memcpy(xx0n,x,elsize);
    for(i=0;i<n;i++)
      memcpy(x+i*elsize,x+(i+1)*elsize,elsize);
    memcpy(x+(n-1)*elsize,xx0n,elsize);
  } else {
    memcpy(xx0n,x+(n-1)*elsize,elsize);
    for(i=(n-1);i>=1;i--)
      memcpy(x+i*elsize,x+(i-1)*elsize,elsize);
    memcpy(x,xx0n,elsize);
  }
  if(xx0n)free(xx0n);
  return;
}
/************************************************************************/
void polar2cart(INT dim, REAL *px, REAL *cx)
{
  // polar is r, theta1,...theta[n-1]; cart is x[0]...x[n-1] px are
  // n-coordinates in polar coord system converts polar coordnates to
  // cartesian in d dimensions.  origin is set at 0,0,0 so if it is
  // different, translation needs to be done after return from here.
  INT i,j;
  REAL rho = px[0];
  REAL cend=rho; 
  switch(dim){
  /* case 1: */
  /*   cx[0]=px[0]; */
  /*   return; */
  /* case 2: */
  /*   cx[0]=rho*cos(px[1]); */
  /*   cx[1]=rho*sin(px[1]); */
  /*   return; */
  /* case 3: */
  /*   cx[0]=rho*sin(px[1])*cos(px[2]); */
  /*   cx[1]=rho*sin(px[1])*sin(px[2]); */
  /*   cx[2]=rho*cos(px[1]); */
  /*   return; */
  default:
    memset(cx,0,dim*sizeof(REAL));
    for(i=0;i<(dim-1);i++){
      cx[i]=rho*cos(px[i+1]);
      for(j=0;j<i;j++){
	cx[i]*=sin(px[j+1]);
	//      fprintf(stdout,"\niiiii=%d,jjjjj=%d",i,j+1);
      }
      //    print_full_mat(1,dim,cx,"c1");
      cend*=sin(px[i+1]);    
    }
    cx[dim-1]=cend;
    // the conversion above puts cx[n-1] first, so put it back at the end.  
    coord_perm(1,dim,cx,sizeof(REAL));
    //        print_full_mat(1,dim,cx,"cx");
    //2d    fprintf(stdout,"\nxc=[%e  %e]",px[0]*cos(px[1]),px[0]*sin(px[1]));
    //    fprintf(stdout,"\nxc=[%e  %e  %e]",			\
    //	    px[0]*sin(px[1])*cos(px[2]),			\
    // px[0]*sin(px[1])*sin(px[2]),				\
    //	    px[0]*cos(px[1])					\
    //	    );
    return;
  }
  return;
}
/************************************************************************/
INT cart2polar(INT dim, REAL *c,REAL *p)
{
  INT i,dimm1=dim-1;
  REAL rl,r;
  // first put c[n-1] first to agree with the polar ordering;
  coord_perm(0,dim,c,sizeof(REAL));
  r=0.;
  for(i=0;i<dim;i++){
    r+=(c[i]*c[i]);
    p[i]=0e0;
  }  
  if(fabs(r)<1e-14){
    for(i=1;i<dim;i++) p[i]=-1e20;
    return 1;
  }
  r=sqrt(r);
  rl=r;
  INT flag=1;
  for(i=1;i<dim;i++){
    p[i]=acos(c[i-1]/rl);
    if(fabs(sin(p[i]))<1e-14){
      flag=0;
      break;
    }    
    rl/=sin(p[i]);
  }
  if(flag) p[dimm1]=atan2(c[dimm1],c[dimm1-1]);
  p[0]=r;
  // permute c back;
  coord_perm(1,dim,c,sizeof(REAL));
  return 0;
}
/************************************************************************/
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void map2mac(scomplex *sc,cube2simp *c2s, input_grid *g)
{
  /* 
     maps a uniform grid from the n-dimensional cube to a hexagonal 
     macroelement from an initial grid of macroelements 
     given by its coordinates xmac[1:nvcube*dim]
     xmac[nvcube][dim]
  */
  INT i,j,k1,k2,k1c,kf,dim=sc->n;
    INT k2c;
  INT ksys;
  REAL *xmac=g->xv;  
  REAL *xhat = (REAL *)calloc(dim,sizeof(REAL));
  REAL *xemac=(REAL *)calloc(c2s->ne*dim,sizeof(REAL));
  // convert midpoints from polar to cartesian.
  //  print_full_mat(c2s->nvcube,c2s->n,g->xv,"P");
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->csysv[k1]];
    k2c=g->systypes[g->csysv[k2]];
    if(g->csysv[k1]==g->csysv[k2] && k1c==1){
      //use xhat as a temp array:
      xhat[0]=0.5*(xmac[k1*dim]+xmac[k2*dim]);// this is rho
      // take half angles;
      for(j=1;j<dim;j++) {
	xhat[j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
      }
      print_full_mat(1,dim,xhat,"xhat");
      polar2cart(dim,xhat,xemac+(i*dim));
      // translate by adding the origin. 
      ksys=g->csysv[k1];// k1c and k2c should be the same below. 
			     k2c=g->csysv[k2];
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
  // now everything is in cartesian, check if  xe are in the convex hull of xmac;
  //midpoints that are
  // not yet attended to are just averages.
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->csysv[k1]];
    k2c=g->systypes[g->csysv[k2]];
    //check all  mid points in polar and skip them; if necessary invert the sign...
    if(g->csysv[k1]==g->csysv[k2] && k1c==1) {      
      REAL scpr1=0.,scpr2=0.;
      for(j=0;j<dim;j++) {
	scpr1+=(xemac[i*dim+j]-g->ox[ksys*dim+j])*(xmac[k1*dim+j]-g->ox[ksys*dim+j]);
	scpr2+=(xemac[i*dim+j]-g->ox[ksys*dim+j])*(xmac[k2*dim+j]-g->ox[ksys*dim+j]);
      }
      // reverse sign if one of the projections on the two vectors forming the arc is negative.
      if(scpr1<0e0 || scpr2<0e00) {
	//	fprintf(stdout,"\nee(%d,%d:%d)=[%d,%d];",i,1,2,k1,k2);
	//	fprintf(stdout,"\n(sc(%d,%d:%d)=(%.16e,%.16e)",i,1,2,scpr1,scpr2);
	for(j=0;j<dim;j++)
	  xemac[i*dim+j]=2e0*g->ox[ksys*dim+j] - xemac[i*dim+j];
      }
      fprintf(stdout,"\ncart:verts=(%d,%d); coord_sys=(%d,%d)",k1,k2,k1c,k2c);fflush(stdout);
      continue;
    }
    for(j=0;j<dim;j++) {
      xemac[i*dim+j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
    }
  }
  print_full_mat(c2s->nvcube,dim,xmac,"X");
  print_full_mat(c2s->ne,dim,xemac,"XE");
  r2c(c2s->nvcube,dim,sizeof(REAL),xmac); // we need xmac by rows here
  r2c(c2s->ne,dim,sizeof(REAL),xemac); // we need xemac (mid points of
				       // edges) also by rows
  for(kf=0;kf<sc->nv;kf++){
    for(i=0;i<dim;i++)xhat[i]=sc->x[kf*dim+i];
    for(i=0;i<dim;i++){
      sc->x[kf*dim+i]=interp8(c2s,xmac+i*c2s->nvcube,xemac+i*c2s->ne,xhat);
      //sc->x[kf*dim+i]=interp4(c2s,xmac+i*c2s->nvcube,xhat);
    }
    //    for(i=0;i<dim;i++){
      //      sc->x[kf*dim+i]=interp4(c2s,xmac+i*c2s->nvcube,xhat);
    //    }
  }
  //  r2c(dim,c2s->nvcube,sizeof(REAL),xmac); // we need xmac by columns here
  //  r2c(dim,c2s->ne,sizeof(REAL),xemac); // we need xemac by rows agin
  if(xhat) free(xhat);
  if(xemac) free(xemac);
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
