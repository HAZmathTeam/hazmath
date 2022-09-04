/*! \file src/amr/amr_coords.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note: containing all routines for polar to cartesian; cartesian to
 *  polar, etc. in any spatial dimension. It also has the routine
 *  mapping "quadratically" the unit square to a macroelement.  Works in dimension n.
 *
 * \todo test the polar and add  cyllindrical and other more general coords. (Ludmil)
 *
 * \note created 20190327 (ltz)
 * \note modified 20190727 (ltz)
 */
#include "hazmath.h"
/**********************************************************************/
/*decaring pi for this scope here*/
static long double pi=M_PI;
/**********************************************************************/
/*!
 * \fn REAL deg2rad(REAL alpha_deg)
 *
 * \param alpha_deg   angle in degrees
 *
 * \return alpha in radians
 *
 */
REAL deg2rad(REAL alpha_deg)
{
  // degrees to radians;
  return (REAL ) ((long double )alpha_deg)*(pi/180e00);
}
/**********************************************************************/
/*!
 * \fn REAL zero_twopi(REAL alpha)
 *
 * \brief map alpha in radians from [-pi,pi) to [0,2*pi]
 *
 * \param alpha   angle in [-pi,pi) in radians
 *
 * \return alpha in [0,2*pi)
 *
 */
REAL zero_twopi(REAL alpha)
{
  long double angle=(long double )alpha;
  angle=atan2l(sinl(angle),cosl(angle));
  if(fabsl(angle)<1e-15) angle=0e0;
  if(angle<0e0)
    return (REAL )(angle+pi+pi);
  return (REAL )angle;
}
/**********************************************************************/
/*!
 * \fn REAL zero_twopi_deg(REAL alpha_deg)
 *
 * \brief map alpha (in degrees) from [-pi,pi) to [0,2*pi) (in radians)
 *
 * \param alpha   angle in [-pi,pi) in degrees
 *
 * \return alpha in [0,2*pi) (in radians)
 *
 */
REAL zero_twopi_deg(REAL alpha_deg)
{
  return zero_twopi(deg2rad(alpha_deg));
}
/**********************************************************************/
/*!
 *\fn static void coord_perm(SHORT type, INT n,void *x, size_t elsize)
 *
 * \brief permutes coordinate arrays (or any array whose elements
 * are of size): if type=1: x[0],x[1]...x[n-1]-->
 * x[1]...x[n-1],x[0]; if type=0(inverse): y[0],y[1]...y[n-1]-->
 * y[n-1],y[0]...y[n-2]
 *
 */
static void coord_perm(const SHORT type, INT n,void *x, size_t elsize)
{
  INT i,n1=n-1;
  void *xx0n=(void *)calloc(1,elsize*sizeof(void));  
  if(type){
    memcpy(xx0n,x,elsize);
    for(i=0;i<n1;i++){
      memcpy(x+i*elsize,x+(i+1)*elsize,elsize);
    }
    memcpy(x+n1*elsize,xx0n,elsize);
  } else {
    memcpy(xx0n,x+n1*elsize,elsize);
    for(i=n1;i>=1;i--)
      memcpy(x+i*elsize,x+(i-1)*elsize,elsize);
    memcpy(x,xx0n,elsize);
  }
  if(xx0n)free(xx0n);
  return;
}
/**********************************************************************/
/*!
 *\fn polar2cart(INT dim, REAL *px, REAL *cx)
 *
 * \brief Converts POLAR (px) TO CARTESIAN(cx) coordinates in any spatial
 *        dimension; Here on input rho=px[0] and px[1],...,px[n-1] are
 *        the angles with px[1] in [0,2*pi) and px[k] in [0,pi],
 *        k=2:(n-1); On output the cartesian coords are x[0]...x[n-1].
 *        The origin is set at (0,...,0).
 *
 *        \param
 *
 * \return
 *
 * \note
 */
/************************************************************************/
void polar2cart(INT dim, REAL *px, REAL *cx)
{
  INT i,j;
  REAL rho = px[0];
  REAL cend=rho;
  /* /\* 1D *\/ */
  /*   cx[0]=px[0]; */
  /*   return; */
  /* /\* 2D *\/ */
  /*   cx[0]=rho*cos(px[1]); */
  /*   cx[1]=rho*sin(px[1]); */
  /*   return; */
  /*   /\*  3D*\/ */
  /*   cx[0]=rho*sin(px[1])*cos(px[2]); */
  /*   cx[1]=rho*sin(px[1])*sin(px[2]); */
  /*   cx[2]=rho*cos(px[1]); */
  /*   return; */
  memset(cx,0,dim*sizeof(REAL));
  for(i=0;i<(dim-1);i++){
    cx[i]=rho*cos(px[i+1]);
    for(j=0;j<i;j++){
      cx[i]*=sin(px[j+1]);
    }
    cend*=sin(px[i+1]);
  }
  cx[dim-1]=cend;
  // the conversion above puts cx[n-1] first, so put it back at the end.
  coord_perm(1,dim,cx,sizeof(REAL));
  return;
}
/*Z!ZZ!Z!ZZ!Z!!Z!Z!Z*/
/**********************************************************************/
/*!
 *\fn INT cart2polar(INT dim, REAL *c,REAL *p)
 *
 * \brief Converts CARTESIAN(c) TO POLAR(p) coordinates in any spatial
 *        dimension; Here, on OUTPUT rho=px[0] and px[1],...,px[n-1]
 *        are the angles with px[1] in [0,2*pi) and px[k] in [0,pi],
 *        k=2:(n-1);
 *
 *
 * \param
 *
 * \return
 *
 * \note
 */
INT cart2polar(INT dim, REAL *c,REAL *p)
{
  /*
   */
  INT i,dimm1=dim-1;
  REAL rl,r;
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
    // BUG  rl/=sin(p[i]);
    rl*=sin(p[i]);
  }
  if(flag) p[dimm1]=atan2(c[dimm1],c[dimm1-1]);
  p[0]=r;
  // permute c back;
  coord_perm(1,dim,c,sizeof(REAL));
  return 0;
}
/**********************************************************************/
/*!
 *\fn void map2mac(scomplex *sc,cube2simp *c2s, input_grid *g)
 *
 *  \brief Maps n-dimensional cube to a macroelement specified by "g".
 *
 * \param
 *
 * \return
 *
 * \todo make all added (refinement) vertices on a boundary that is in polar coordinates to be
 * in the same polar coordinate system  (ltz).
 */
void map2mac(scomplex *sc,cube2simp *c2s, input_grid *g)
{
  /*
  */
  INT i,j,k1,k2,k1c,kf,dim=sc->n;
  //INT k2c;
  INT ksys;
  REAL *xmac=g->xv;
  REAL *xhat = (REAL *)calloc(dim,sizeof(REAL));
  REAL *c1 = (REAL *)calloc(dim,sizeof(REAL));
  REAL *c2 = (REAL *)calloc(dim,sizeof(REAL));
  REAL *xemac=(REAL *)calloc(c2s->ne*dim,sizeof(REAL));
  REAL rho;
  //  input_grid_print(g);
  // convert midpoints from polar to cartesian.
  //  print_full_mat(c2s->nvcube,c2s->n,g->xv,"P");
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->csysv[k1]];
    //k2c=g->systypes[g->csysv[k2]];
    if(g->csysv[k1]==g->csysv[k2] && k1c==1){
      //use xhat as a temp array:
      rho=0.5*(xmac[k1*dim]+xmac[k2*dim]);// this is the rho we will use
      /*
	 To find the mid point in polar: convert the vertices k1 and
	 k2 to cartesian, take the middle point (in cartesian) and
	 then use the angles defined for the middle point and average
	 radius rho
      */
      polar2cart(dim,xmac+k1*dim,c1); polar2cart(dim,xmac+k2*dim,c2);
      for(j=0;j<dim;j++)
	c1[j]=0.5*(c1[j]+c2[j]);
      cart2polar(dim,c1,xhat);
      //      print_full_mat(1,dim,xhat,"xhat");
      xhat[0]=rho;
      polar2cart(dim,xhat,xemac+(i*dim));
      // translate by adding the origin.
      ksys=g->csysv[k1];// k1c and k2c should be the same below.
      //k2c=g->csysv[k2];
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
  //midpoints that are not yet attended to are just averages.
  for(i=0;i<c2s->ne;i++){
    k1=c2s->edges[2*i];
    k2=c2s->edges[2*i+1];
    k1c=g->systypes[g->csysv[k1]];
    //k2c=g->systypes[g->csysv[k2]];
    //skip all  mid points in polar
    if(g->csysv[k1]==g->csysv[k2] && k1c==1) continue;
    for(j=0;j<dim;j++) {
      xemac[i*dim+j]=0.5*(xmac[k1*dim+j]+xmac[k2*dim+j]);
    }
  }
  //  print_full_mat(c2s->nvcube,dim,xmac,"X");
  //  print_full_mat(c2s->ne,dim,xemac,"XE");
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
  free(xhat);
  free(xemac);
  free(c1);
  free(c2);
  return;
}
/*EOF*/
